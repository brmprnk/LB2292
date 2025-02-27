import xml.etree.ElementTree as ET
import os
import pandas as pd
import pprint
import os


files = []
dirs = os.listdir("project2/data/raw/clinical")
for dir in dirs:
    for file in os.listdir("project2/data/raw/clinical/" + dir):
        if file.endswith(".xml"):
            # append the full path, starting at project2/
            files.append("project2/data/raw/clinical/" + dir + "/" + file)


def rec_parse(root):
    def rec_call(el):
        if len(el) == 0:
            return el.text
        else:
            return {el.tag.split('}')[-1]: rec_call(el) for el in el}
    return root.tag.split('}')[-1], rec_call(root)


def parse_new_tumor_events(el):
    """
    Parses the new tumor events section of the clinical data.
    Since this is stored differently for melanoma, metastasis and neo-plasm events, we have to combine them here a bit.

    Parameters:
    el (Element): The element to parse.

    Returns:
    dict: The parsed dictionary.
    
    """
    parsed = {}

    for nte in el:
        res = parse_flatten(nte, prefix="recurrence")
        if res == {}:
            continue
    
        # These are the keys that are always in there:
        recurrence_cols = [
            'recurrence_new_neoplasm_event_type',
            'recurrence_days_to_new_tumor_event_after_initial_treatment',
            'recurrence_additional_radiation_therapy',
            'recurrence_additional_pharmaceutical_therapy',
            'recurrence_additional_surgery_locoregional_procedure',
            'recurrence_additional_surgery_metastatic_procedure',

            # The following were not included, but could be if we really want to (not super useful I think)
            # 'recurrence_nte_pathologic_tumor_length',
            # 'recurrence_nte_pathologic_tumor_width',
            # 'recurrence_nte_pathologic_tumor_depth',
            # 'recurrence_nte_radiologic_tumor_length',
            # 'recurrence_nte_radiologic_tumor_width',
            # 'recurrence_nte_radiologic_tumor_depth',
        ]


        parsed["recurrence_new_tumor_anatomic_site"] = None
        if "recurrence_new_neoplasm_event_occurrence_anatomic_site" in res.keys() and res["recurrence_new_neoplasm_event_occurrence_anatomic_site"] is not None:
            parsed["recurrence_new_tumor_anatomic_site"] = res["recurrence_new_neoplasm_event_occurrence_anatomic_site"]
            if res["recurrence_new_neoplasm_event_occurrence_anatomic_site"].lower() == "other, specify" and "recurrence_new_neoplasm_occurrence_anatomic_site_text" in res.keys():
                parsed["recurrence_new_tumor_anatomic_site"] = res["recurrence_new_neoplasm_occurrence_anatomic_site_text"]
            # print("tumor site:", parsed["recurrence_new_tumor_anatomic_site"])
        elif "recurrence_new_tumor_metastasis_anatomic_site" in res.keys() and res["recurrence_new_tumor_metastasis_anatomic_site"] is not None:
            parsed["recurrence_new_tumor_anatomic_site"] = res["recurrence_new_tumor_metastasis_anatomic_site"]
            if res["recurrence_new_tumor_metastasis_anatomic_site"].lower() == "other, specify" and "recurrence_new_tumor_metastasis_anatomic_site_other_text" in res.keys():
                parsed["recurrence_new_tumor_anatomic_site"] = res["recurrence_new_tumor_metastasis_anatomic_site_other_text"]
            # print("metastasis site:", parsed["recurrence_new_tumor_anatomic_site"])
        elif "recurrence_new_primary_melanoma_anatomic_site" in res.keys() and res["recurrence_new_primary_melanoma_anatomic_site"] is not None:
            parsed["recurrence_new_tumor_anatomic_site"] = res["recurrence_new_primary_melanoma_anatomic_site"]

        for k in recurrence_cols:
            if k in res.keys():
                parsed[k] = res[k]
            else:
                parsed[k] = None

        # rename recurrence_new_neoplasm_event_type to recurrence_new_tumor_event_type
        if "recurrence_new_neoplasm_event_type" in parsed.keys():
            parsed["recurrence_new_tumor_event_type"] = parsed["recurrence_new_neoplasm_event_type"]
            del parsed["recurrence_new_neoplasm_event_type"]
        
        # if recurrence_new_non_melanoma_event_histologic_type_text is filled, and not "No", we put it in the event type column
        if "recurrence_new_non_melanoma_event_histologic_type_text" in parsed.keys() and parsed["recurrence_new_non_melanoma_event_histologic_type_text"] is not None and parsed["recurrence_new_non_melanoma_event_histologic_type_text"].lower() != "no":
            parsed["recurrence_new_tumor_event_type"] = parsed["recurrence_new_non_melanoma_event_histologic_type_text"]

        # capitalize only the first letter of each word in "recurrence_new_tumor_anatomic_site"
        if "recurrence_new_tumor_anatomic_site" in parsed.keys() and parsed["recurrence_new_tumor_anatomic_site"] is not None:
            parsed["recurrence_new_tumor_anatomic_site"] = parsed["recurrence_new_tumor_anatomic_site"].title()
            # TODO: maybe I want to do this for all non-numerical columns actually?

        # if there is anything else in res that is not in kvs.keys(), print those entries
        # for k in res.keys():
        #     if k not in recurrence_cols and res[k] is not None:
        #         print("Unused key:", k, res[k])
    return parsed


all_responses = []


def parse_drugs(el):

    # Relevant fields:
    keep_cols = [
        # "drug_name",  # (too many different ones and combinations, so we don't use this)
        "therapy_type",
        "therapy_type_notes",
        "therapy_ongoing",  
        "measure_of_response",
    ]

    parsed = {}
    parsed["drug_received_treatment"] = (len(el) >= 1)  # True if any drug is found

    received_treatments = []
    responses = []
    therapy_ongoing = None

    # if len(el) > 10:
    for drug in el:
        drug = parse_flatten(drug, prefix="")

        # Filter only relevant columns
        drug = {k: v for k, v in drug.items() if k in keep_cols}
        # Put all missing to None
        for col in keep_cols:
            if col not in drug.keys():
                drug[col] = None

        # Rescue some of the "other" values in therapy_type, otherwise set to "Other"
        if drug["therapy_type"] is not None and drug["therapy_type"].lower() == "other, specify in notes":
            if drug["therapy_type_notes"] in ["Ancillary", "Ancillary agent", "Ancillary Agent"]:
                drug["therapy_type"] = "Ancillary"
            elif drug["therapy_type_notes"] in ["Vaccine", "vaccine"]:
                drug["therapy_type"] = "Vaccine"
            else:
                drug["therapy_type"] = "Other"

        if drug["therapy_type"] is None:
            drug["therapy_type"] = "None"  # TODO: check when this occurs. It shouldn't?

        therapy_type = drug["therapy_type"].lower().replace(" ", "_")

        received_treatments.append(drug["therapy_type"])
        responses.append(drug["measure_of_response"])
        all_responses.append(drug["measure_of_response"])

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


def parse_clinical_list(tag: str, column: ET.Element) -> tuple:
    """
    Parse function for if the clinical data field is a list-type structure.
    """
    values = [v.text for v in column]
    if values == [None]:
        return tag, pd.NA
    else:
        return tag, ", ".join(values)
    

def parse_leaf_values(el, join=False):
    def rec_call(x):
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


def rec_find(el, tag, last=False):
    # last: if True, return the last value found without warning
    value = None
    for child in el:
        if child.tag.split('}')[-1] == tag:
            if value is not None and value != child.text:
                print(f"Warning, multiple values found for tag={tag}, value={value}, new_value={child.text}")
            value = child.text
        _, child = rec_find(child, tag)
        if child is not None:
            if value is not None and value != child:
                print(f"Warning, multiple values found for tag={tag}, value={value}, new_value={child}")
            value = child
    return tag, value


def parse_relative_cancer_history_list(el):

    relatives = []
    for relative in el:
        d = {}
        for c in el[0]:
            subtag = c.tag.split("}")[-1]
            subvalue = c.text
            d[subtag] = subvalue

        type_col = "cancer_diagnosis_cancer_type_icd9_text_name"
        if not type_col in d.keys():
            type_col = "relative_family_cancer_hx_text"
        if not type_col in d.keys():
            type_col = "family_history_cancer_type_other"
        if not type_col in d.keys():
            type_col = "family_cancer_type_txt"
        try:
            comb = f'{d["family_member_relationship_type"]}:{d[type_col]}'
        except:
            print(f"Unable to parse family history: {d}")
        if comb is not None and comb != "None:None" and comb not in relatives:
            relatives.append(comb)

    if len(relatives) > 1:
        print(f"Warning, multiple relatives found, not supported: {relatives}")

    if relatives == []:
        value = None
    else:
        value = ", ".join(relatives)

    return "blood_relative_cancer_history_list", value


def parse_flatten(el, prefix: str = "") -> dict:
    """
    Flattens the contents of this elements and returns it as a dictionary. Prefix is prepended to all keys of the leaf nodes.
    
    Parameters:
    el (Element): The element to parse
    prefix (str): The prefix to prepend to the keys of the leaf nodes.

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
    return parsed



def parse_dev(el):
    print(el.tag)
    pprint.pp(rec_parse(el))
    return el.tag, None


def parse_file(FILE):
    parsed = {}

    admin, patient = ET.parse(FILE).getroot()

    _, disease_code = rec_find(admin, "disease_code")
    parsed["cancer_type"] = disease_code

    for child in patient:

        key = child.tag.split('}')[-1]
        value = None

        if len(child) == 0:  # Easy cases, just key-value pairs.
            key, value = key, child.text

        elif key in ["race_list",
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
            _, parsed["ablation_performed_indicator"] = rec_find(child, "ablation_performed_indicator", last=True)
            types = rec_find(child, "ablation_type")
            if types is not None and type(types) == list and len(types) > 1:
                types = ", ".join(types)
            parsed["ablation_type"] = types
            _, parsed["emolization_performed_indicator"] = rec_find(child, "emolization_performed_indicator", last=True)
            types = rec_find(child, "emolization_type")
            if types is not None and type(types) == list and len(types) > 1:
                types = ", ".join(types)
            parsed["emolization_type"] = types
            continue
            
        elif key == "neoplasm_dimension":
            _, parsed["neoplasm_length"] = rec_find(child, "neoplasm_length", last=True)
            _, parsed["neoplasm_width"] = rec_find(child, "neoplasm_width", last=True)
            _, parsed["neoplasm_depth"] = rec_find(child, "neoplasm_depth", last=True)
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
            # example: {'molecular_analysis_abnormality_testing_result_values': {'molecular_analysis_abnormality_testing_result': 'NPMc Negative',
            #    'molecular_analysis_abnormality_testing_percentage_value': None}})
            # the percentage is always None, so we only use testing_result
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
            key, value = parse_relative_cancer_history_list(child)
            if value is None:
                continue
            relation, cancer = value.split(":")
            parsed["blood_relative_cancer_history_relation"] = relation
            parsed["blood_relative_cancer_history_cancertype"] = cancer
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

        elif key == "radiations":
            continue  # Not used, because we already have the "radiation_therapy" column, which is just binary, but sufficient.

        elif key == "treatment":
            continue  # No additional information in here

        elif key == "history_of_disease":
            for k, v in parse_flatten(child, prefix="history").items():
                parsed[k] = v
            continue

        elif key == "tests_performed":
            continue  # NOTE: Cool information, but very specific, so not likely to be useful for students.

        elif key == "stage_event":
            # Useful information about the cancer stage (not comparable between cancer types though)
            for k, v in parse_flatten(child, prefix="stage").items():
                parsed[k] = v
            continue

        elif key == "follow_ups":
            """
            This is actually available for 86% of the patients, whereas primary_therapy_outcome is missing for 80%
            Can I use this to fill in some of the missing values? (i.e. Does the followup contain outcome information if the primary therapy outcome is missing?)
            """


            for follow_up in child:
                flat =  parse_flatten(follow_up, prefix="follow_up")

                if "follow_up_primary_therapy_outcome_success" in flat.keys():
                    parsed["follow_up_primary_therapy_outcome_success"] = flat["follow_up_primary_therapy_outcome_success"]
                    # NOTE: this is only 59% missing in the followup data, instead of 80% in the primary therapy outcome data

                if "follow_up_vital_status" in flat.keys():
                    parsed["follow_up_vital_status"] = flat["follow_up_vital_status"]
                    # NOTE: this went from 75% alive to %61 alive in the followup data
                    # Oh, but this is because the followup data is more recent (i.e. people may have died of other causes in between), so it's not that useful.

                if "follow_up_followup_treatment_success" in flat.keys():
                    parsed["follow_up_followup_treatment_success"] = flat["follow_up_followup_treatment_success"]
                    # NOTE: this is 75% missing

                # Hmm, I'm not sure how to use this exactly.
                # How about...
                # - Use follow-up data to fill out the primary therapy outcome if it's missing
                # But what if they disagree? Do I take the follow-up because it's more recent? 
                # We have follow_up_followup_treatment_success, follow_up_primary_therapy_outcome_success and primary_therapy_outcome_success
                #     and they are not the same in 510 rows
                # And what do I do with multuple follow-ups???

                # primary_pathology_first_treatment_success --> this one is useless I think
                
                # pprint.pp(flat)

            # print("===================\n\n")

            key, value = "follow_ups", "Available"  # TODO: do we need this? Does this contain extra survival data?

        elif key == "anatomic_neoplasm_subdivisions":  # Not important, hard to parse
            continue
        elif key == "first_nonlymph_node_metastasis_anatomic_sites":
            continue  # Not needed
        elif key == "additional_studies":
            continue  # Doesn't seem useful.
        elif key == "patient_history_immune_system_and_related_disorders_names":
            continue  # Only 3 values filled (out of which two are "other"), so we skip it.
        elif key == "genotyping_results_gene_mutation_not_reported_reasons":
            continue  # Not needed
        elif key == "history_non_muscle_invasive_blca":
            continue  # No clue what this is, so I'm skipping it for now.
        elif key == "FISH_test_component_results":
            # example: {'FISH_test_component_result': {'FISH_test_component': 'PML-RAR',
                                #  'FISH_test_component_percentage_value': '95'}})
            # however, it may also have multiple values. NOTE: skipped for now.
            continue
        else:
            print(f"Warning, unimplemented key: {key}")

        if value is not None:
            parsed[key] = value

    return parsed


if __name__ == "__main__":

    parsed_files = []

    from collections import defaultdict

    cancertype_counts = defaultdict(int)

    for file in files:
        parsed = parse_file(file)
        # cancertype_counts[parsed["cancer_type"]] += 1
        parsed_files.append(parsed)

    df = pd.DataFrame(parsed_files)
    # df = df.filter(regex='^rad').dropna(how='all')  # Only for debugging the recurrence data

    # TODO: filter a lot of columns (we have like 750 now, and most of them are mostly empty)
    # Filter such that:
    # For each column there should be at least one cancer type for which the column is 50% filled, otherwise drop it.

    print(df.head())

    # to csv
    df.to_csv("project2/data/clinical_parsed_recurrence.csv", index=False)