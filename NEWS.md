# EcoTools 0.1.0

## New features
* Added pre-fit fixed-effect design-matrix validation.
* Added sparse-response screening for count and binary models.
* Added analysis fingerprints for cache validation.
* Added repeated DHARMa diagnostic-stability assessment.
* Added standardized model-status and inference-validity fields.
* Added support for validated estimated marginal means.

## Tweedie simulation safety
* Added a pre-simulation safety guard for fitted Tweedie models.
* DHARMa simulation is skipped when the maximum fitted Tweedie power is at or above 1.99 by default.
* The guard applies to supported `glmmTMB` Tweedie models and `mgcv` Tweedie GAMs.
* For `mgcv::twlss()` models, the guard evaluates the maximum fitted power across observations before allocating or generating the simulation matrix.
* Repeated DHARMa stability assessment is skipped for models that trigger the guard.
* Guarded model bundles record the fitted power range, configured threshold, machine-readable reason, warning message, and DHARMa status.
* Models that trigger the guard are retained for inspection but are not marked valid for inference.
The threshold can be configured with `tweedie_simulation_max_power`.

## Bug fixes
* Prevented Tweedie simulations near the Gamma boundary from starting before the simulation-safety check.
* Initialized Tweedie guard and DHARMa status metadata on all fitting paths, including non-Tweedie models and failed fits.
* Included the Tweedie simulation threshold in cache fingerprints.
* Improved direct-model support in `DHARMa_diag()`.
* Improved serialization of model objects fitted with function-local data.