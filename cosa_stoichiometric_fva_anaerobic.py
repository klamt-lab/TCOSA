from cosa_stoichiometric_fva import cosa_single_swap_test
import ray
ray.init(log_to_driver=False)

for base_nadx_scenario in ("WILDTYPE", "FLEXIBLE", "SINGLE_COFACTOR"):
    cosa_single_swap_test(
        anaerobic=True,
        reac_id="TEST_0_321",
        mu=0.321,
        base_nadx_scenario=base_nadx_scenario
    )
