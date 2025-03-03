import os

import pandas as pd
import xml.etree.ElementTree as ET

from tqdm import tqdm

from process_utils import find_data_files


def parse_recursive(root: ET.Element) -> dict:
    """
    Recusively parses the XML tree and returns it as a dictionary.
    Note: if there is nested structure, this will also return nested dictionaries.
        For a flat structure, use parse_flatten instead.

    Parameters:
    root (Element): The root element to parse.

    Returns:
    dict: The parsed dictionary.
    """
    def rec_call(el):
        if len(el) == 0:
            return el.text
        else:
            return { child.tag.split('}')[-1]: rec_call(child) for child in el }
        
    return root.tag.split('}')[-1], rec_call(root)


def parse_new_tumor_events(el: ET.Element) -> dict:
    """
    Parses the new tumor events section of the clinical data.
    Since this is stored differently for melanoma, metastasis and neo-plasm events, we have to combine them here a bit.

    Parameters:
    el (Element): The element to parse.

    Returns:
    dict: The parsed dictionary.
    """
    nte_types = []
    nte_sites = []

    for nte in el:
        res = parse_flatten(nte, prefix="nte")
        if res == {}:
            continue

        for type_col in [  # Get the type from one of these columns, and put it in the event type column
            "nte_new_neoplasm_event_type",
            "nte_new_non_melanoma_event_histologic_type_text",
            "nte_new_tumor_anatomic_site",
        ]:
            if type_col in res.keys() and res[type_col] is not None:
                nte_types.append(res[type_col])

        for site_col in [  # Get the site from one of these columns, and put it in the event site column
            "nte_new_neoplasm_event_occurrence_anatomic_site",
            "nte_new_neoplasm_occurrence_anatomic_site_text",
            "nte_new_tumor_metastasis_anatomic_site",
            "nte_new_tumor_metastasis_anatomic_site_other_text",
            "nte_new_primary_melanoma_anatomic_site",
        ]:
            if site_col in res.keys() and res[site_col] is not None and res[site_col].lower() != "other, specify":
                nte_sites.append(res[site_col])

        # Not using the following, because they are hard to combine if there are multiple events
        # 'nte_days_to_new_tumor_event_after_initial_treatment',
        # 'nte_additional_radiation_therapy',
        # 'nte_additional_pharmaceutical_therapy',
        # 'nte_additional_surgery_locoregional_procedure',
        # 'nte_additional_surgery_metastatic_procedure',

        # The following were not included, but could be if we really want to (not super useful I think)
        # 'nte_nte_pathologic_tumor_length',
        # 'nte_nte_pathologic_tumor_width',
        # 'nte_nte_pathologic_tumor_depth',
        # 'nte_nte_radiologic_tumor_length',
        # 'nte_nte_radiologic_tumor_width',
        # 'nte_nte_radiologic_tumor_depth',

    return {
        "new_tumor_event_type": nte_types[0] if len(set(nte_types)) == 1 else "Multiple",
        "new_tumor_event_site": nte_sites[0] if len(set(nte_sites)) == 1 else "Multiple",
    }


def parse_drugs(el: ET.Element) -> dict:
    """
    Parse the drugs section of the clinical data.
    Most of the data is very specific, so we only parse the therapy_type, therapy_ongoing and measure_of_response fields.
    The therapy_type_notes field is used to rescue some of the "other" values in therapy_type.
    We add columns for each therapy type (boolean for whether the patient recieved this treatment),
    and for each therapy type we add a column for the response.

    Parameters:
    el (Element): The element to parse.

    Returns:
    dict: The parsed dictionary.
    """
    parsed = { "drug_received_treatment": len(el) >= 1 }  # True if any drug is found
    received_treatments = []
    responses = []
    therapy_ongoing = None

    for drug in el:
        # Only keep the relevant columns, and set missing to None
        drug = parse_flatten(drug, prefix="", filter=["therapy_type", "therapy_type_notes", "therapy_ongoing", "measure_of_response"])

        # Rescue some of the "other" values in therapy_type, otherwise set to "Other"
        if drug["therapy_type"] is not None and drug["therapy_type"].lower() == "other, specify in notes":
            if drug["therapy_type_notes"] in ["Ancillary", "Ancillary agent", "Ancillary Agent"]:
                drug["therapy_type"] = "Ancillary"
            elif drug["therapy_type_notes"] in ["Vaccine", "vaccine"]:
                drug["therapy_type"] = "Vaccine"
            else:
                drug["therapy_type"] = "Other"

        if drug["therapy_type"] is None:
            continue  # Skip this entry if the therapy type is missing

        therapy_type = drug["therapy_type"].lower().replace(" ", "_")

        received_treatments.append(drug["therapy_type"])
        responses.append(drug["measure_of_response"])

        # We set "received_treatment_xxx" to True if any drug of that type is found. If multiple exist, we only keep the last one.
        parsed[f"drug_received_treatment_{therapy_type}"] = True

        # But for the response, we want update the response if it's better, such that we store only the best response.
        def better_reponse(new_response, old_response):
            order = [None, 'Clinical Progressive Disease', 'Stable Disease', 'Partial Response', 'Complete Response']
            return order.index(new_response) > order.index(old_response)
        if f"drug_received_treatment_{therapy_type}_response" not in parsed.keys() \
                or better_reponse(drug["measure_of_response"], parsed[f"drug_received_treatment_{therapy_type}_response"]):
            parsed[f"drug_received_treatment_{therapy_type}_response"] = drug["measure_of_response"]

        # Finally, we want to store whether ANY treatment is still ongoing, so we store False when we find False, unless True was already there, and always write True
        if drug["therapy_ongoing"] is not None and drug["therapy_ongoing"].upper() == "YES":
            therapy_ongoing = True
        elif drug["therapy_ongoing"] is not None and drug["therapy_ongoing"].upper() == "NO" and therapy_ongoing is None:
            therapy_ongoing = False

    parsed["drug_therapy_ongoing"] = therapy_ongoing
    return parsed
    

def parse_leaf_values(el: ET.Element, join=False) -> tuple:
    """
    Parses the leaf values of the element, and returns them as a list.
    So this is similar to parse_flatten, but only for the leaf nodes.

    Parameters:
    el (Element): The element to parse.
    join (bool): If True, join the values with a comma.

    Returns:
    tuple: The tag and the value.
    """
    def rec_call(x):  # Recursive call 
        if len(x) == 0:
            return [x.text]
        else:
            ret = []
            for y in x:
                ret.extend(rec_call(y))
        
    key = el.tag.split('}')[-1]
    if len(el) == 0:
        value = el.text
    else:
        values = []
        for child in el:
            rec_res = rec_call(child)
            if rec_res is not None:
                values.extend(rec_res)
        values = [v for v in values if v is not None]

        if values == [None] or values == [] or ", ".join(values) == "None":
            value = None
        elif join:
            value = ", ".join(sorted(values))
        else:
            value = values  # Only then it's returned as a list
    return key, value


def rec_find(el: ET.Element, tag: str) -> tuple:
    """
    Recursively finds the tag in the element, and returns the value.

    Parameters:
    el (Element): The element to parse.
    tag (str): The tag to find.

    Returns:
    tuple: The tag and the value.
    """
    value = None
    for child in el:
        if child.tag.split('}')[-1] == tag:
            if value is not None and value != child.text:
                print(f"Warning, multiple values found for tag={tag}, value={value}, new_value={child.text}")
            value = child.text
        _, child_val = rec_find(child, tag)
        if child_val is not None:
            if value is not None and value != child_val:
                print(f"Warning, multiple values found for tag={tag}, value={value}, new_value={child_val}")
            value = child_val
    return tag, value


def parse_relative_cancer_history_list(el: ET.Element) -> dict:
    """
    Parses the blood_relative_cancer_history_list section of the clinical data.
    This is a list of family members who have had cancer, and the type of cancer they had.
    However, we simplify into "Mutiple" for both relation and type if there are multiple (different) values.
    So we always end up with two columns.

    Parameters:
    el (Element): The element to parse.

    Returns:
    dict: The parsed dictionary.
    """
    relationship_types = []
    cancer_types = []

    for relative in el:
        d = { c.tag.split("}")[-1]: c.text for c in relative } 

        # if all values are None, continue
        if all([v is None for v in d.values()]):
            continue

        # Cancer type may be denoted in several columns, we pick the first one we find (there's usually only one filled)
        possible_type_cols = ["cancer_diagnosis_cancer_type_icd9_text_name",
                              "relative_family_cancer_hx_text",
                              "family_history_cancer_type_other",
                              "family_history_cancer_type",
                              "family_cancer_type_txt"]
        for type_col in possible_type_cols:
            if type_col in d.keys() and d[type_col] is not None:
                break
        if type_col not in d.keys() or d[type_col] is None or "family_member_relationship_type" not in d.keys() or d["family_member_relationship_type"] is None:
            continue  # If one of them is None, skip this entry

        relationship_types.append(d["family_member_relationship_type"].replace("Natural ", ""))  # (remove prefix Natural, same thing)
        cancer_types.append(d[type_col].title())

    if len(relationship_types) == 0:
        return {}

    return {
        "blood_relative_cancer_history_relation": relationship_types[0] if len(relationship_types) == 1 else "Multiple",
        "blood_relative_cancer_history_cancertype": cancer_types[0] if len(cancer_types) == 1 else "Multiple",
    }


def parse_flatten(el, prefix: str = "", filter=None) -> dict:
    """
    Flattens the contents of this elements and returns it as a dictionary. Prefix is prepended to all keys of the leaf nodes.
    If a filter list is provided, only the fields in the filter list are kept, or set to None if they are missing.
    
    Parameters:
    el (Element): The element to parse
    prefix (str): The prefix to prepend to the keys of the leaf nodes.
    filter (list): A list of keys to keep. If a key is not in this list, it is set to None.

    Returns:
    dict: The parsed dictionary.
    """
    if prefix != "" and not prefix.endswith("_"):  # Make sure the prefix ends with an underscore
        prefix = prefix + "_"

    parsed = {}
    for child in el:
        if len(child) == 0:
            parsed[prefix + child.tag.split('}')[-1]] = child.text
        else:  # Recurse
            for k, v in parse_flatten(child, prefix).items():
                parsed[k] = v

    if filter is not None:
        parsed = {k: v for k, v in parsed.items() if k in filter}  # Keep only the fields in the keep list
        for col in filter:  # Put all missing to None
            if col not in parsed.keys():
                parsed[col] = None

    return parsed


def parse_file(file_path: str) -> dict:
    """
    Parses a single clinical file, at the given file path.
    Delegates to the more specific parsing functions, based on what is required for that field.
    This is also where we decide what columns to parse, and which to skip. 
    For more complex fields, an explanation of the approach is given in the corresponding case below.

    Parameters:
    file_path (str): The path to the file to parse.

    Returns:
    dict: The parsed dictionary.
    """
    parsed = {}

    admin, patient = ET.parse(file_path).getroot()

    # From the admin header, we only need the disease code (the cancer type)
    _, disease_code = rec_find(admin, "disease_code")
    parsed["cancer_type"] = disease_code

    for child in patient:
        key = child.tag.split('}')[-1]
        value = None

        if key in [
            "anatomic_neoplasm_subdivisions",
            "first_nonlymph_node_metastasis_anatomic_sites",
            "additional_studies",
            "patient_history_immune_system_and_related_disorders_names",
            "genotyping_results_gene_mutation_not_reported_reasons",
            "history_non_muscle_invasive_blca",
            "FISH_test_component_results",
            "radiations",  # Not used, because we already have the "radiation_therapy" column, which is just binary, but sufficient.
            "treatment",  # No additional information in here
            "tests_performed",  # Not used, because it's very specific, and not likely to be useful for students. (But it contains some cool things)
        ]:
            continue  # These are the columns that are not used.

        elif len(child) == 0:  # Easy cases, just key-value pairs.
            key, value = key, child.text

        elif key in [
            "race_list",
            "viral_hepatitis_serologies",
            "first_degree_relative_history_thyroid_gland_carcinoma_diagnosis_relationship_types",
            "lymph_node_preoperative_assessment_diagnostic_imaging_types",
            "metastatic_neoplasm_confirmed_diagnosis_method_names",
            "sites_of_primary_melanomas",
            "metastatic_site_list",
            "tumor_levels",
            "diagnostic_mri_results",
            "antireflux_treatment_types",
            "cytogenetic_abnormalities",
            "relation_testicular_cancer_list",
            "lymph_node_location_positive_pathology_names",
            "human_papillomavirus_types",
            "fdg_or_ct_pet_performed_outcomes",
            "relative_cancer_types",
            "patient_personal_medical_history_thyroid_gland_disorder_names",
        ]:
            # These fields are technically lists, but it rarely happens that there are multiple values, 
            #   so we just join them with a comma when that happens.
            key, value = parse_leaf_values(child, join=True)
        
        elif key == "loss_expression_of_mismatch_repair_proteins_by_ihc_results":
            # This just has a group around it that we don't need, so we can just flatten it
            key, values = parse_leaf_values(child)
            for v in ['MLH1-Expressed', 'MSH2-Expressed', 'PMS2-Expressed', 'MSH6-Expressed']:
                if values is not None and v in values:
                    parsed[f"{key}_{v}"] = True
            for v in ['MLH1-Not Expressed', 'MSH2-Not Expressed', 'PMS2-Not Expressed', 'MSH6-Not Expressed']:
                if values is not None and v in values:
                    v = v.replace("Not ", "")  # (make it the same column)
                    parsed[f"{key}_{v}"] = False
            continue

        elif key == "history_hepato_carcinoma_risk_factors":
            key, value = parse_leaf_values(child, join=True)
            if value is not None:
                value = value.replace("Other, ", "").replace(", Other", "").replace("Other", "")  # remove the "Other" values
                if value == "":
                    value = None

        elif key == "ablations":
            # NOTE: there is a bit more information here about the procedure, but I think only these four are enough
            # The type fields can have several values, so I joined them with a comma for now.
            _, parsed["ablation_performed_indicator"] = rec_find(child, "ablation_performed_indicator")
            types = rec_find(child, "ablation_type")
            if types is not None and type(types) == list and len(types) > 1:
                types = ", ".join(types)
            parsed["ablation_type"] = types
            _, parsed["emolization_performed_indicator"] = rec_find(child, "emolization_performed_indicator")
            types = rec_find(child, "emolization_type")
            if types is not None and type(types) == list and len(types) > 1:
                types = ", ".join(types)
            parsed["emolization_type"] = types
            continue
            
        elif key == "neoplasm_dimension":
            _, parsed["neoplasm_length"] = rec_find(child, "neoplasm_length")
            _, parsed["neoplasm_width"] = rec_find(child, "neoplasm_width")
            _, parsed["neoplasm_depth"] = rec_find(child, "neoplasm_depth")
            continue        

        elif key == "prior_systemic_therapy_types":
            key, value = rec_find(child, "prior_systemic_therapy_type")

        elif key == "diagnostic_ct_abd_pelvis_results":
            key, value = rec_find(child, "diagnostic_ct_abd_pelvis_result")
        
        elif key == "immunophenotype_cytochemistry_testing_results":
            parsed["immunophenotype_cytochemistry_testing_performed"] = True
            key, values = parse_leaf_values(child)
            if values is None:
                continue
            # Example of values: [['MPX Positive', None], ['NSE Negative', None], ['TDT Not Tested', None], ['CD3 Not Tested', None], ['CD5 Not Tested', None], ['CD7 Not Tested', None], ['CD10 Not Tested', None], ['CD14 Not Tested', None], ['CD15 Not Tested', None], ['CD19 Not Tested', None], ['CD20 Not Tested', None], ['CD23 Not Tested', None], ['CD25 Not Tested', None], ['CD33 Positive', None], ['CD34 Not Tested', None], ['CD45 Not Tested', None], ['CD56 Not Tested', None], ['CD117 Positive', None], ['HLA-DR Positive', None], ['Other CD Not Tested', None]]
            # We only want the first value of each pair, and we want to remove the "Not Tested" values.
            values = [v[0] for v in values if v[0]]
            # Now all end in "Positive", "Negative", or "Not Tested", so we make a column for each marker, and fill this as the value.
            for v in values:
                if v.endswith("Positive"):
                    parsed[f"immunophenotype_cytochemistry_testing_{v.replace(' Positive', '')}"] = "Positive"
                elif v.endswith("Negative"):
                    parsed[f"immunophenotype_cytochemistry_testing_{v.replace(' Negative', '')}"] = "Negative"
                elif v.endswith("Not Tested"):
                    parsed[f"immunophenotype_cytochemistry_testing_{v.replace(' Not Tested', '')}"] = "Not Tested"
            continue

        elif key == "molecular_analysis_abnormality_testing_results":
            parsed["molecular_analysis_abnormality_testing_performed"] = True
            key, values = parse_leaf_values(child)
            if values is None:
                continue
            # Example of values: [['MPX Positive', None], ['NSE Negative', None], ['TDT Not Tested', None], ['CD3 Not Tested', None], ['CD5 Not Tested', None], ['CD7 Not Tested', None], ['CD10 Not Tested', None], ['CD14 Not Tested', None], ['CD15 Not Tested', None], ['CD19 Not Tested', None], ['CD20 Not Tested', None], ['CD23 Not Tested', None], ['CD25 Not Tested', None], ['CD33 Positive', None], ['CD34 Not Tested', None], ['CD45 Not Tested', None], ['CD56 Not Tested', None], ['CD117 Positive', None], ['HLA-DR Positive', None], ['Other CD Not Tested', None]]
            # We only want the first value of each pair, and we want to remove the "Not Tested" values.
            values = [v[0] for v in values if v[0] is not None]
            # Now all end in "Positive", "Negative", so we make a column for each marker, and fill this as the value.
            for v in values:
                if v.endswith("Positive"):
                    parsed[f"molecular_analysis_abnormality_testing_{v.replace(' Positive', '')}"] = "Positive"
                elif v.endswith("Negative"):
                    parsed[f"molecular_analysis_abnormality_testing_{v.replace(' Negative', '')}"] = "Negative"
            continue

        elif key in ["diagnostic_ct_result_outcomes", "diagnostic_mri_result_outcomes"]:
            key, values = parse_leaf_values(child)
            if values is None:
                continue            
            # All values end in Present or Absent, so we can make a column for each.
            for v in values:
                if v.endswith("Present"):
                    parsed[f"{key}_{v.replace(' Present', '')}"] = True
                elif v.endswith("Absent"):
                    parsed[f"{key}_{v.replace(' Absent', '')}"] = False
            continue

        elif key == "postoperative_tx_list":
            key, value = rec_find(child, "postoperative_tx")

        elif key == "blood_relative_cancer_history_list":
            for k, v in parse_relative_cancer_history_list(child).items():
                parsed[k] = v
            continue

        elif key == "primary_pathology":
            for k, v in parse_flatten(child, prefix="primary_pathology").items():
                parsed[k] = v
            continue

        elif key == "new_tumor_events":
            for k, v in parse_new_tumor_events(child).items():
                parsed[k] = v
            continue

        elif key == "drugs":
            for k, v in parse_drugs(child).items():
                parsed[k] = v
            continue

        elif key == "history_of_disease":
            for k, v in parse_flatten(child, prefix="history").items():
                parsed[k] = v
            continue

        elif key == "stage_event":
            # Useful information about the cancer stage (not comparable between cancer types though)
            for k, v in parse_flatten(child, prefix="stage").items():
                parsed[k] = v
            continue

        elif key == "follow_ups":
            # There is a lot more information here, but we only use the last available therapy outcome and vital status.
            # In post-processing, we combine the follow_up data with the "primary_therapy_outcome" column
            #    to get the "final_therapy_outcome", containing the most recent data.
            for follow_up in child:
                flat = {k: v for k, v in parse_flatten(follow_up, prefix="follow_up").items() if v is not None}
                if "follow_up_primary_therapy_outcome_success" in flat.keys():
                    parsed["follow_up_primary_therapy_outcome_success"] = flat["follow_up_primary_therapy_outcome_success"]
                if "follow_up_vital_status" in flat.keys():
                    parsed["follow_up_vital_status"] = flat["follow_up_vital_status"]
                if "follow_up_followup_treatment_success" in flat.keys():
                    parsed["follow_up_treatment_success"] = flat["follow_up_followup_treatment_success"]
                continue

        else:
            print(f"Warning, unimplemented key: {key}")

        if value is not None:
            parsed[key] = value

    return parsed


def parse_clinical(path: str, all_columns: bool = False) -> pd.DataFrame:
    """
    Parses all clinical files into a single DataFrame.
    Then we filter the columns to get rid of a lot of the mostly-empty ones.

    We filter such that:
    - For each column there should be at least one cancer type for which the column is 50% filled, otherwise drop it.
    - Also remove columns that are missing for more than 95% of the samples (regardless of cancer type)
    - Combine 

    Parameters:
    all_columns (bool): If True, keep all columns, otherwise only keep the most important ones

    Returns:
    pd.DataFrame: The parsed DataFrame.
    """
    df = pd.DataFrame([parse_file(file) for file in tqdm(find_data_files(path, ext=".xml"), desc="Parsing clinical data")])

    # Add new column that combines the treatment outcomes such that we keep the most recent one.
    df['final_treatment_outcome'] = df['follow_up_treatment_success'].combine_first(
        df['follow_up_primary_therapy_outcome_success']).combine_first(df['primary_therapy_outcome_success'])
    
    if all_columns:
        return df
    
    # Drop these, because they were combined into the final_treatment_outcome column
    df = df.drop(columns=['follow_up_treatment_success', 'follow_up_primary_therapy_outcome_success', 'primary_therapy_outcome_success'])

    # Only keep columns if there is at least one cancer type for which the column is filled for at least 50% of the samples
    keep_columns = set()
    for cancer_type in sorted(df["cancer_type"].unique()):
        # add columns that are filled for at least 50% for this cancer type to the keep_columns list
        df_cancer_type = df[df["cancer_type"] == cancer_type]
        missing_values = df_cancer_type.isnull().mean(axis=0)
        keep_columns.update(missing_values[missing_values < 0.5].index.tolist())
    keep_columns = list(keep_columns)

    # Also remove columns that are missing for more than 95% of the samples
    missing_values = df.isnull().mean(axis=0)
    keep_columns = [col for col in keep_columns if missing_values[col] < 0.95]

    # Also include all the columns I put too much effort in while parsing
    keep_columns += [col for col in df.columns if "drug_" in col]
    keep_columns += [col for col in df.columns if "blood_relative_" in col]

    # Remove duplicates and sort alphabetically
    keep_columns = sorted(list(set(keep_columns)))
    df = df[keep_columns].copy()

    # rename bcr_patient_barcode to patient_id, and set it as the index
    df["patient_id"] = df["bcr_patient_barcode"]
    df.set_index("bcr_patient_barcode", inplace=True)

    return df


if __name__ == "__main__":

    OUT_PATH = "project2/data/processed/clinical_full.csv"
    df = parse_clinical(path="project2/data/raw/clinical")
    df.to_csv(OUT_PATH, index=False)

    print(df)
    print(f"Saved clinical data to {OUT_PATH}")
