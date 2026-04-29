using MLJ
using MLJLinearModels
using EvoTrees
import MLJDecisionTreeInterface

"""
    superlearner(; metalearner=nothing, binary=false)

Return an MLJ stacked ensemble (SuperLearner).
`binary=true` for propensity score (classifier), `binary=false` for outcome regression.
"""
function superlearner(; metalearner=nothing, binary=false)
    if binary
        RF   = MLJDecisionTreeInterface.RandomForestClassifier
        meta = metalearner === nothing ? LogisticClassifier() : metalearner
        return Stack(
            metalearner = meta,
            resampling  = CV(nfolds=5),
            glm  = LogisticClassifier(),
            rf   = RF(n_trees=500),
            evo  = EvoTreeClassifier(nrounds=100, max_depth=5),
            mean = ConstantClassifier()
        )
    else
        RF   = MLJDecisionTreeInterface.RandomForestRegressor
        meta = metalearner === nothing ? LinearRegressor() : metalearner
        return Stack(
            metalearner = meta,
            resampling  = CV(nfolds=5),
            glm  = LinearRegressor(),
            rf   = RF(n_trees=500),
            evo  = EvoTreeRegressor(nrounds=100, max_depth=5),
            mean = ConstantRegressor()
        )
    end
end

function _mlj_fit_classifier(model, X::DataFrame, y)
    mach = machine(model, X, categorical(y .== 1))
    MLJ.fit!(mach, verbosity=0)
    mach
end

function _mlj_predict_prob1(mach, X::DataFrame)
    preds = MLJ.predict(mach, X)
    clip(Float64[pdf(p, true) for p in preds], 1e-8, 1 - 1e-8)
end

function _mlj_fit_regressor(model, X::DataFrame, y)
    mach = machine(model, X, Float64.(y))
    MLJ.fit!(mach, verbosity=0)
    mach
end

function _mlj_predict_reg(mach, X::DataFrame)
    Float64.(MLJ.predict(mach, X))
end
