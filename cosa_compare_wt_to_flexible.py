import cobra
from helper import json_zip_load, json_load, get_files
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids

dG0_values = json_load("cosa/dG0_values.json")
cobra_model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")

def print_reacs(reacs):
    for reac in reacs:
        cobra_reac = cobra_model.reactions.get_by_id(reac)
        try:
            dG0 = round(dG0_values[reac]["dG0"], 2)
        except KeyError:
            dG0 = float("NaN")
        print(f"* {reac};  ΔG0={dG0} kJ/mol; {cobra_reac.reaction}")

tcosa_reactions = get_all_tcosa_reaction_ids(cobra_model)

for aerobicity in ("aerobic",):# "anerobic"):
    print(f"={aerobicity}=")
    for concentrations in ("STANDARDCONC",): # "VIVOCONC"):
        for target in ("OPTSUBMDF",):  # "OPTMDF"
            nums_common = []
            nums_uncommon = []
            print(f"=={target}==")
            wt_data = json_zip_load(
                f"cosa/results_{aerobicity}/runs/{target}_{concentrations}_WILDTYPE.json"
            )["0,568"]
            current_sampling = 0
            for file in get_files(f"cosa/results_{aerobicity}/flux_samplings"):
                print("~~~", current_sampling, "~~~")
                current_sampling += 1
                fx_data = json_load(
                    f"cosa/results_{aerobicity}/flux_samplings/{file}"
                )

                essential_reacs = {'KARA2_FWD_ORIGINAL_NADP_TCOSA', 'HISTD_ORIGINAL_NAD_TCOSA', 'TYRL_ORIGINAL_NADP_TCOSA', 'EAR40y_ORIGINAL_NADP_TCOSA', 'DPR_ORIGINAL_NADP_TCOSA', 'EGMEACPR_ORIGINAL_NADP_TCOSA', 'EPMEACPR_ORIGINAL_NADP_TCOSA', 'AMPMS2_ORIGINAL_NAD_TCOSA', 'AGPR_REV_ORIGINAL_NADP_TCOSA', 'HSDy_REV_ORIGINAL_NADP_TCOSA', 'KARA1_REV_ORIGINAL_NADP_TCOSA', '3OAR80_FWD_ORIGINAL_NADP_TCOSA', 'EAR60y_ORIGINAL_NADP_TCOSA', 'SULR_ORIGINAL_NADP_TCOSA', '3OAR140_FWD_ORIGINAL_NADP_TCOSA', 'DHDPRy_ORIGINAL_NADP_TCOSA', '3OAR40_FWD_ORIGINAL_NADP_TCOSA', 'DHFR_FWD_ORIGINAL_NADP_TCOSA', 'P5CR_ORIGINAL_NADP_TCOSA', 'PPND_ORIGINAL_NAD_TCOSA', 'NADS1_ORIGINAL_NAD_TCOSA', 'PDX5PS_ORIGINAL_NAD_TCOSA', 'ICDHyr_FWD_ORIGINAL_NADP_TCOSA', 'GLYCL_ORIGINAL_NAD_TCOSA', 'UAPGR_ORIGINAL_NADP_TCOSA', 'GLUTRR_ORIGINAL_NADP_TCOSA', 'SHK3Dr_FWD_ORIGINAL_NADP_TCOSA', 'MTHFR2_ORIGINAL_NAD_TCOSA', 'G5SD_ORIGINAL_NADP_TCOSA', '3OAR100_FWD_ORIGINAL_NADP_TCOSA', 'OGMEACPR_ORIGINAL_NADP_TCOSA', 'THZPSN3_ORIGINAL_NADP_TCOSA', '3OAR60_FWD_ORIGINAL_NADP_TCOSA', 'IPMD_ORIGINAL_NAD_TCOSA', 'EAR80y_ORIGINAL_NADP_TCOSA', 'OPMEACPR_ORIGINAL_NADP_TCOSA', 'APRAUR_ORIGINAL_NADP_TCOSA', 'DXPRIi_ORIGINAL_NADP_TCOSA'}
                found_essential_reacs = []
                num_common = 0
                num_uncommon = 0
                # for growth_rate in wt_data.keys():
                wt_values = wt_data["values"]
                fx_values = fx_data["values"]
                uncommon_reacs = []
                for tcosa_reaction in tcosa_reactions:
                    # if ("_ORIGINAL_" in tcosa_reaction) or (tcosa_reaction in ["NADK", "NADTRHD", "NADPPPS", "THD2pp"]):
                    #     pass
                    # else:
                    #     continue

                    if tcosa_reaction not in essential_reacs:
                        continue

                    z_var = f"{tcosa_reaction}"
                    try:
                        wt_values[z_var]
                    except KeyError:
                        continue
                    on_in_wt = wt_values[z_var] > 1e-7
                    on_in_fx = fx_values[z_var] > 1e-7


                    off_in_wt = abs(wt_values[z_var]) < 1e-7
                    off_in_fx = abs(fx_values[z_var]) < 1e-7

                    if (off_in_fx and off_in_wt):
                        continue

                    if (on_in_fx and on_in_wt):
                        num_common += 1
                    else:
                        num_uncommon += 1
                        uncommon_reacs.append(tcosa_reaction)

                    found_essential_reacs.append(tcosa_reaction)

                nums_common.append(num_common)
                nums_uncommon.append(num_uncommon)
                print(f"Common: {num_common}; uncommon: {num_uncommon}")
                not_found = []
                for essential_reac in essential_reacs:
                    if essential_reac not in found_essential_reacs:
                        not_found.append(essential_reac)
                # print(" ->Uncommon:")
                # print_reacs(uncommon_reacs)
                # print(" ->Not found essential:", not_found)

            print(sum(nums_common) / len(nums_common))
            print(sum(nums_uncommon) / len(nums_uncommon))
