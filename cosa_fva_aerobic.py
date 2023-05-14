from cosa_fva import cosa_single_swap_test
import ray
ray.init(log_to_driver=False)

for base_nadx_scenario in ("WILDTYPE", "FLEXIBLE"):
    cosa_single_swap_test(
        anaerobic=False,
        reac_id="TEST_0_868",
        mu=0.868,
        base_nadx_scenario=base_nadx_scenario
    )
    cosa_single_swap_test(
        anaerobic=False,
        reac_id="TEST_0_818",
        mu=0.818,
        base_nadx_scenario=base_nadx_scenario
    )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_768",
    #     mu=0.768,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_718",
    #     mu=0.718,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    cosa_single_swap_test(
        anaerobic=False,
        reac_id="TEST_0_668",
        mu=0.668,
        base_nadx_scenario=base_nadx_scenario
    )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_618",
    #     mu=0.618,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_568",
    #     mu=0.568,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_518",
    #     mu=0.518,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_468",
    #     mu=0.468,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    cosa_single_swap_test(
        anaerobic=False,
        reac_id="TEST_0_418",
        mu=0.418,
        base_nadx_scenario=base_nadx_scenario
    )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_368",
    #     mu=0.418,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_418",
    #     mu=0.418,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_368",
    #     mu=0.368,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_318",
    #     mu=0.318,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_268",
    #     mu=0.268,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    cosa_single_swap_test(
        anaerobic=False,
        reac_id="TEST_0_218",
        mu=0.218,
        base_nadx_scenario=base_nadx_scenario
    )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_168",
    #     mu=0.168,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_118",
    #     mu=0.118,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_068",
    #     mu=0.068,
    #     base_nadx_scenario=base_nadx_scenario
    # )
    cosa_single_swap_test(
        anaerobic=False,
        reac_id="TEST_0_05",
        mu=0.05,
        base_nadx_scenario=base_nadx_scenario
    )
    # cosa_single_swap_test(
    #     anaerobic=False,
    #     reac_id="TEST_0_03",
    #     mu=0.03,
    #     base_nadx_scenario=base_nadx_scenario
    # )
