classdef BayesoptOptions_custom
    %

    %   Copyright 2016-2020 The MathWorks, Inc.

    properties
        % Properties that can be set using bayesopt name/value pairs. Any
        % defaults here define the name/value pair defaults.
        
        % Controlling the search:
        AcquisitionFunctionName = 'expectedimprovementpersecondplus';
        IsObjectiveDeterministic = false;
        ExplorationRatio = 0.5;
        
        % Parallel
        UseParallel = false;
        ParallelMethod = 'clippedmodelprediction';
        MinWorkerUtilization = [];                                          % Default depends on pool size.
        
        % Starting and Stopping:
        NumSeedPoints = 4;
        MaxObjectiveEvaluations = [];                                       % Default: 30 (or gridsize if grid)
        MaxTime = BayesoptOptions_custom.DefaultMaxTime();
        
        % Fitting
        GPActiveSetSize = 300;
        
        % Constrained optimization:
        XConstraintFcn = [];
        ConditionalVariableFcn = [];
        NumCoupledConstraints = 0;
        CoupledConstraintTolerances = [];
        AreCoupledConstraintsDeterministic = [];
       
        % Reporting, Saving, and Interrupting execution:
        Verbose = 1;
        OutputFcn = {};
        SaveVariableName = 'BayesoptResults';
        SaveFileName = 'BayesoptResults.mat';

        % Plotting:
        PlotFcn = {@plotObjectiveModel, @plotMinObjective};
        
        % Initialization from prior function evaluations:
        InitialX = [];
        InitialObjective = [];
        InitialConstraintViolations = [];
        InitialErrorValues = [];
        InitialObjectiveEvaluationTimes = [];
        InitialIterationTimes = [];
        InitialUserData = {};
    end
    
    properties(Dependent)
        VariableDescriptions;
    end
    
    properties(Hidden)
        % Problem definition
        ObjectiveFcn;
        NumGridDivisions = 10;
        ObjectiveNargout = [];
        
        % Main algorithm:
        AlwaysReportObjectiveErrors = false;
        DefaultMinWorkerUtilizationFraction = 0.8;
        
        % Constrained Optimization:
        ErrorConstraintTol = .01;
        XConstraintSatisfiabilitySamples = 10000;
        
        % Plotting and Display:
        NumPlotGrid = [];
        ShowEnvelopes = [];
        DisplayHeaderInterval = 20;
        IsClassregRegressionFunction = false;
        
        % Function modeling:
        FitModels = true;
        Sigma0Divisor = 5;
        SigmaFTol = 1e-3;
        MaxExploitIterations = 5;
        CatKernelScale = 0.2;
        
        % Optimizing the acquisition function:
        NumRestarts = 10;
        NumRestartCandidates = 1e4;
        VerboseRestarts = false;
        MaxIterPerRestart = 10;
        RelTol = 1e-3;   
        
        % The type of surrogate model used. Supports 'SingleTreeBagger' and
        % 'MultiTreeBagger' for Random Forest Regression and
        % 'GaussianProcess' for Gaussian Process Regression.
        % MultiTreeBagger is currently hard-coded to split on an a variable
        % called 'algorithm'.
        ModelType = 'GaussianProcess';
        
        % Num of initial points to evaluate for each algorithm when
        % ModelType is 'MultiTreeBagger'.
        NumInitialTrialsForEachAlgo = 10;
        
        % boolean set to true if bayesopt is called from 'fitcauto'.
        FitcautoMultipleLearners = false;
        
        % boolean set to true if bayesopt is called from 'fitcauto'.
        FitcAutoSingleLearner = false;
    end
    
    properties(Hidden, Dependent)
        VarSpec;
    end        
    
    properties(Access=protected)
        PrivVarSpec;
        PrivVariableDescriptions;
    end
    
    methods
        % Constructor
        function this = BayesoptOptions_custom(ObjectiveFcn, VariableDescriptions, NVPs)
            this.ObjectiveFcn = ObjectiveFcn;
            this.VariableDescriptions = VariableDescriptions;
            this = fillFromNVPs(this, NVPs);
            this = validateAndFillDefaults(this);
        end
        
        function this = fillFromNVPs(this, RemainingNVPs)
            % Parse all user-visible public properties of BayesoptOptions_custom,
            % plus these hidden properties: FitModels, NumGridDivisions,
            % AlwaysReportObjectiveErrors, IsClassregRegressionFunction.
            Names = [properties(this); {'FitModels'; 'NumGridDivisions'; 'AlwaysReportObjectiveErrors'; 'IsClassregRegressionFunction'; 'ModelType';'FitcautoMultipleLearners';'FitcAutoSingleLearner'}];
            [Vals{1:numel(Names)}, setflags, RemainingNVPs] = internal.stats.parseArgs(Names, repmat({[]},numel(Names),1), RemainingNVPs{:});
            for i = 1:numel(Names)
                if setflags.(Names{i})
                    this.(Names{i}) = Vals{i};
                end
            end
            if ~isempty(RemainingNVPs)
                bayesoptim.err('UnknownBayesoptArg', RemainingNVPs{1});
            end
        end
        
        function this = validateAndFillDefaults(this)
            % Validate user-supplied options and fill defaults. Order is
            % important.
            this = checkAndFillParallelComputing(this);
            checkSaveVariableName(this);
            checkSaveFileName(this);
            checkVerbose(this);
            checkModelType(this);
            checkFitcautoMultipleLearners(this);
            checkNumCoupledConstraints(this);
            this = checkAndFillNumGridDivisions(this);
            checkFitModels(this);
            checkAlwaysReportObjectiveErrors(this);
            checkIsClassregRegressionFunction(this);
            this = checkObjectiveFcn(this);
            this = checkXConstraintFcn(this);
            this = checkConditionalVariableFcn(this);
            this = checkAndFillSearchParams(this);
            this = checkAndFillStoppingCriteria(this);
            checkGPActiveSetSize(this);
            checkNumSeedPoints(this);
            checkNumInitialTrialsForEachAlgo(this);
            this = checkAndFillConstraints(this);
            this = fillPlotParams(this);
            this = checkAndFillPlotFcn(this);
            checkOutputFcn(this);
            this = checkAndFillInitializationData(this);
        end
        
        % Setters and getters
        function this = set.VariableDescriptions(this, VariableDescriptions)
            this.PrivVariableDescriptions = VariableDescriptions;
            this.PrivVarSpec = iCreateVarSpecFromBayesoptVars(VariableDescriptions);
        end
        
        function VariableDescriptions = get.VariableDescriptions(this)
            VariableDescriptions = this.PrivVariableDescriptions;
        end
        
        function VarSpec = get.VarSpec(this)
            VarSpec = this.PrivVarSpec;
        end
    end
    
    methods(Access=protected)
        function checkSaveVariableName(this)
            if ~isvarname(this.SaveVariableName)
                bayesoptim.err('SaveVariableName');
            end
        end
        
        function checkSaveFileName(this)
            if ~(ischar(this.SaveFileName) && isvector(this.SaveFileName))
                bayesoptim.err('SaveFileName');
            end
        end
        
        function checkVerbose(this)
            if ~bayesoptim.isLowerBoundedIntScalar(this.Verbose, 0)
                bayesoptim.err('Verbose');
            end
        end
        
        function checkModelType(this)
            % Check the given option is one of 'GaussianProcess',
            % 'SingleTreeBagger' or 'MultiTreeBagger'
            if ~ismember(this.ModelType, {'GaussianProcess', 'MultiGaussianProcess', 'MultiTreeBagger', 'SingleTreeBagger', 'RandomSearch'})
                bayesoptim.err('ModelType',this.ModelType);
            end
        end
        
        function checkFitcautoMultipleLearners(this)
            % Check the given option is a boolean
            validateattributes(this.FitcautoMultipleLearners,{'logical'},{'scalar'})
        end
        
        function checkFitcAutoSingleLearner(this)
            % Check the given option is a boolean
            validateattributes(this.FitcAutoSingleLearner,{'logical'},{'scalar'})
        end
        
        function checkNumCoupledConstraints(this)
            if ~bayesoptim.isLowerBoundedIntScalar(this.NumCoupledConstraints, 0)
                bayesoptim.err('NumCoupledConstraintsType');
            end
        end
        
        function this = checkObjectiveFcn(this)
            % Check datatype
            if this.UseParallel
                if ~(isa(this.ObjectiveFcn, 'function_handle') || ...
                     isa(this.ObjectiveFcn, 'parallel.pool.Constant'))
                    bayesoptim.err('ObjectiveFcnTypeParallel');
                end
            else
                if ~isa(this.ObjectiveFcn, 'function_handle')
                    bayesoptim.err('ObjectiveFcnTypeSerial');
                end
            end
            % Get nargout now if possible
            if isa(this.ObjectiveFcn, 'function_handle')
                Nargout = getObjFcnNargout(this.ObjectiveFcn, this.NumCoupledConstraints);
                if Nargout > 0
                    this.ObjectiveNargout = Nargout;
                end
            end
        end
        
        function this = checkXConstraintFcn(this)
            % The output of the function gets checked every call.
            if ~isempty(this.XConstraintFcn) && ~isa(this.XConstraintFcn, 'function_handle')
                bayesoptim.err('XConstraintFcnType');
            end
        end
        
        function this = checkConditionalVariableFcn(this)
            % The output of the function gets checked every call.
            if ~isempty(this.ConditionalVariableFcn) && ~isa(this.ConditionalVariableFcn, 'function_handle')
                bayesoptim.err('ConditionalVariableFcnType');
            end
        end
        
        % Controlling the search:
        function this = checkAndFillSearchParams(this)
            if isempty(this.AcquisitionFunctionName) || ~ischar(this.AcquisitionFunctionName)
                bayesoptim.err('AFNameType');
            end
            RepairedValue = bayesoptim.parseArgValue(this.AcquisitionFunctionName, ...
                {'expectedimprovement', 'expectedimprovementpersecond',...
                'lowerconfidencebound', 'probabilityofimprovement',...
                'expectedimprovementplus', 'expectedimprovementpersecondplus',...
                'grid', 'random'});
            if isempty(RepairedValue)
                bayesoptim.err('AFName', this.AcquisitionFunctionName);
            else
                this.AcquisitionFunctionName = RepairedValue;
            end
            if ~bayesoptim.isLogicalScalar(this.IsObjectiveDeterministic)
                bayesoptim.err('IsObjectiveDeterministic');
            end
            if ~bayesoptim.isNonnegativeRealScalar(this.ExplorationRatio)
                bayesoptim.err('ExplorationRatio');
            end
        end
        
        % Parallel computing:
        function this = checkAndFillParallelComputing(this)
            % UseParallel
            if ~bayesoptim.isLogicalScalar(this.UseParallel)
                bayesoptim.err('UseParallel');
            end
            if this.UseParallel && ~hasParallel()
                bayesoptim.err('NoParallelComputingToolbox');
            end
            % ParallelMethod
            if isempty(this.ParallelMethod) || ~ischar(this.ParallelMethod)
                bayesoptim.err('ParallelMethodType');
            end
            RepairedValue = bayesoptim.parseArgValue(this.ParallelMethod, ...
                {'clippedmodelprediction', 'modelprediction', ...
                'minobserved', 'maxobserved'});
            if isempty(RepairedValue)
                bayesoptim.err('ParallelMethod', this.ParallelMethod);
            else
                this.ParallelMethod = RepairedValue;
            end
            % MinWorkerUtilization
            if ~isempty(this.MinWorkerUtilization) && ~bayesoptim.isLowerBoundedIntScalar(this.MinWorkerUtilization, 1)
                bayesoptim.err('MinWorkerUtilization');
            end
        end
        
        % Stopping criteria:
        function this = checkAndFillStoppingCriteria(this)
            if isempty(this.MaxObjectiveEvaluations)
                if isequal(this.AcquisitionFunctionName, 'grid')
                    this.MaxObjectiveEvaluations = prod(this.NumGridDivisions);
                else
                    this.MaxObjectiveEvaluations = BayesoptOptions_custom.DefaultMaxObjectiveEvaluations();
                end
            end
%             if ~bayesoptim.isLowerBoundedIntScalar(this.MaxObjectiveEvaluations,1)
%                 bayesoptim.err('MaxObjectiveEvaluations');
%             end
            if ~bayesoptim.isNonnegativeRealScalar(this.MaxTime)
                bayesoptim.err('MaxTime');
            end
        end
        
        function checkGPActiveSetSize(this)
            if ~bayesoptim.isLowerBoundedIntScalar(this.GPActiveSetSize, 1)
                bayesoptim.err('GPActiveSetSize');
            end
        end
        
        % Constrained optimization:
        function this = checkAndFillConstraints(this)
            if ~isempty(this.XConstraintFcn)
                if ~isa(this.XConstraintFcn, 'function_handle')
                    bayesoptim.err('XConstraintFcnType');
                elseif nargin(this.XConstraintFcn) == 0 || nargout(this.XConstraintFcn) == 0
                    bayesoptim.err('XConstraintFcnNarg');
                end
            end
            if ~isempty(this.ConditionalVariableFcn)
                if ~isa(this.ConditionalVariableFcn, 'function_handle')
                    bayesoptim.err('ConditionalVariableFcnType');
                elseif nargin(this.ConditionalVariableFcn) == 0 || nargout(this.ConditionalVariableFcn) == 0
                    bayesoptim.err('ConditionalVariableFcnNarg');
                end
            end
            if isempty(this.CoupledConstraintTolerances)
                if this.NumCoupledConstraints > 0
                    % fill
                    this.CoupledConstraintTolerances = .01*ones(1, this.NumCoupledConstraints);
                end
            else
                % check type and size
                if ~(isnumeric(this.CoupledConstraintTolerances) && ...
                     isreal(this.CoupledConstraintTolerances) && ...
                     isvector(this.CoupledConstraintTolerances) && ...
                     all(this.CoupledConstraintTolerances >= 0) && ...
                     numel(this.CoupledConstraintTolerances) == this.NumCoupledConstraints)
                    bayesoptim.err('CoupledConstraintTolerances');
                end
            end
            if isempty(this.AreCoupledConstraintsDeterministic)
                if this.NumCoupledConstraints > 0
                    % fill
                    this.AreCoupledConstraintsDeterministic = true(1, this.NumCoupledConstraints);
                end
            else
                % check type and size
                if ~(islogical(this.AreCoupledConstraintsDeterministic) && ...
                        isvector(this.AreCoupledConstraintsDeterministic) && ...
                        numel(this.AreCoupledConstraintsDeterministic) == this.NumCoupledConstraints)
                    bayesoptim.err('AreCoupledConstraintsDeterministic');
                end
            end
        end
        
        % Plotting:
        function this = fillPlotParams(this)
            D = this.VarSpec.NumVars;
            if D==1
                this.NumPlotGrid = 1000;
            elseif D==2
                this.NumPlotGrid = [50, 50];
            end
            this.ShowEnvelopes = D==1;
        end
        
        function this = checkAndFillPlotFcn(this)
            if isempty(this.PlotFcn) || isa(this.PlotFcn, 'function_handle')
                return;
            elseif iscell(this.PlotFcn) && ...
                    all(cellfun(@(x)isa(x, 'function_handle'), this.PlotFcn))
                return;
            elseif isequal(this.PlotFcn, 'all')
                this.PlotFcn = {...
                    @plotObjectiveModel,...
                    @plotAcquisitionFunction,...
                    @plotObjectiveEvaluationTimeModel,...
                    @plotConstraintModels,...
                    @plotObjective,...
                    @plotObjectiveEvaluationTime,...
                    @plotMinObjective,...
                    @plotElapsedTime};
                return;
            else
                bayesoptim.err('PlotFcn');
            end
        end
        
        function checkOutputFcn(this)
            if isempty(this.OutputFcn) || isa(this.OutputFcn, 'function_handle')
                return;
            elseif iscell(this.OutputFcn) && ...
                    all(cellfun(@(x)isa(x, 'function_handle'), this.OutputFcn))
                return;
            else
                bayesoptim.err('OutputFcn');
            end
        end
        
        % Initialization from prior function evaluations:
        function this = checkAndFillInitializationData(this)
            AllNonXEmpty = isempty(this.InitialObjective) && ...
                isempty(this.InitialConstraintViolations) && ...
                isempty(this.InitialErrorValues) && ...
                isempty(this.InitialObjectiveEvaluationTimes) && ...
                isempty(this.InitialIterationTimes) && ...
                isempty(this.InitialUserData);
            if isempty(this.InitialX)
                if AllNonXEmpty
                    % Everything is empty.
                    return;
                else
                    % X is empty but something else isn't.
                    bayesoptim.err('InitialXEmpty');
                end
            else
                % X is nonempty. Check it.
                this = checkInitialX(this);
                if AllNonXEmpty
                    % Only X is provided.
                    return;
                else
                    % X is nonempty and so are some of the others. Check
                    % nonempties and fill empties with NaNs.
                    N = size(this.InitialX,1);
                    this = checkInitialF(this, N);
                    this = checkInitialConstraintViolations(this, N);
                    this = checkInitialErrorValues(this, N);
                    this = checkInitialObjectiveEvaluationTimes(this, N);
                    this = checkInitialIterationTimes(this, N);
                    this = reconcileInitialIterationAndObjectiveTimes(this);
                    this = checkInitialUserData(this, N);
                end
            end
        end
        
        function this = checkInitialX(this)
            % Check type
            if ~istable(this.InitialX)
                bayesoptim.err('InitialXType');
            end
            % Check width
            if size(this.InitialX, 2) ~= this.VarSpec.NumVars
                bayesoptim.err('InitialXWidth', size(this.InitialX, 2), this.VarSpec.NumVars);
            end
            % Check each var
            for i = 1:this.VarSpec.NumVars
                X = this.InitialX{:,i};
                if this.VarSpec.isCat(i)
                    if ~iscategorical(X)
                        bayesoptim.err('InitialXCatType', i);
                    end
                    % Make sure passed cats are known
                    if ~all(ismember(X, this.VarSpec.Categories{i}))
                        bayesoptim.err('InitialXCats', i);
                    end
                    % Convert InitialX categories to VarSpec categories so
                    % their double representations match.
                    this.InitialX.(i) = categorical(X, categories(this.VarSpec.Categories{i}));
                else
                    % This var is numeric
                    if ~bayesoptim.isAllFiniteReal(X)
                        bayesoptim.err('InitialXNumericType', i);
                    end
                    if ~this.VarSpec.isReal(i)
                        % Verify ints were passed
                        if ~all(bayesoptim.isInteger(X))
                            bayesoptim.err('InitialXInts', i);
                        end
                    end
                    % Verify passed log vars are positive
                    if this.VarSpec.isLog(i) && ~all(X > 0)
                        bayesoptim.err('InitialXLogs', i);
                    end
                    % All checks passed. Convert to double
                    this.InitialX.(i) = double(this.InitialX.(i));
                end
            end
        end
        
        function this = checkInitialF(this, N)
            % F must be empty or have the same height as X and be real.
            if isempty(this.InitialObjective)
                this.InitialObjective = NaN(N,1);
            else
                if size(this.InitialObjective,1)~=N
                    bayesoptim.err('InitialObjectiveHeight');
                elseif ~bayesoptim.isAllFiniteRealOrNaN(this.InitialObjective)
                    bayesoptim.err('InitialObjectiveType');
                end
                % All checks passed. Convert to double
                this.InitialObjective = double(this.InitialObjective);
            end
        end
        
        function this = checkInitialConstraintViolations(this, XHeight)
            % ConstraintViolations must be empty or have the same height as X and be
            % real.
            if isempty(this.InitialConstraintViolations)
                this.InitialConstraintViolations = NaN(XHeight, this.NumCoupledConstraints);
            else
                if size(this.InitialConstraintViolations,1) ~= XHeight
                    bayesoptim.err('InitialConstraintViolationsHeight',...
                        size(this.InitialConstraintViolations,1), XHeight);
                elseif size(this.InitialConstraintViolations,2) ~= this.NumCoupledConstraints
                    bayesoptim.err('InitialConstraintViolationsWidth',...
                        size(this.InitialConstraintViolations,2), this.NumCoupledConstraints);
                elseif ~bayesoptim.isAllFiniteRealOrNaN(this.InitialConstraintViolations)
                    bayesoptim.err('InitialConstraintViolationsType');
                end
                % All checks passed. Convert to double
                this.InitialConstraintViolations = double(this.InitialConstraintViolations);
            end
        end
        
        function this = checkInitialObjectiveEvaluationTimes(this, N)
            % InitialObjectiveEvaluationTimes must be empty or have the
            % same height as X, be real and positive or NaN.
            T = this.InitialObjectiveEvaluationTimes;
            if isempty(T)
                this.InitialObjectiveEvaluationTimes = NaN(N,1);
            else
                if size(T,1)~=N
                    bayesoptim.err('InitialObjectiveEvaluationTimesSize');
                elseif ~all(isnumeric(T) && isreal(T))
                    bayesoptim.err('InitialObjectiveEvaluationTimesType');
                elseif ~all(isnan(T) | ~isinf(T) & T>0)
                    bayesoptim.err('InitialObjectiveEvaluationTimesType');
                end
                % All checks passed. Convert to double
                this.InitialObjectiveEvaluationTimes = double(this.InitialObjectiveEvaluationTimes);
            end
        end
        
        function this = checkInitialIterationTimes(this, N)
            % InitialIterationTimes must be empty or have the same height
            % as X, be real and positive or NaN.
            T = this.InitialIterationTimes;
            if isempty(T)
                this.InitialIterationTimes = NaN(N,1);
            else
                if size(T,1)~=N
                    bayesoptim.err('InitialIterationTimesSize');
                elseif ~all(isnumeric(T) && isreal(T))
                    bayesoptim.err('InitialIterationTimesType');
                elseif ~all(isnan(T) | ~isinf(T) & T>0)
                    bayesoptim.err('InitialIterationTimesType');
                end
                % All checks passed. Convert to double
                this.InitialIterationTimes = double(this.InitialIterationTimes);
            end
        end
        
        function this = reconcileInitialIterationAndObjectiveTimes(this)
            % Overwrite iterationtime nans with feval times where they
            % exist.
            IsnanIter = isnan(this.InitialIterationTimes);
            IsnanEval = isnan(this.InitialObjectiveEvaluationTimes);
            this.InitialIterationTimes(IsnanIter & ~IsnanEval) = ...
                this.InitialObjectiveEvaluationTimes(IsnanIter & ~IsnanEval);
        end
        
        function this = checkInitialUserData(this, N)
            % InitialUserData must be empty or be a cell vector the same
            % height as X.
            T = this.InitialUserData;
            if isempty(T)
                this.InitialUserData = cell(N,1);
            else
                if ~iscell(T)
                    bayesoptim.err('InitialUserDataType');
                elseif ~(isvector(T) || length(T)==N)
                    bayesoptim.err('InitialUserDataSize');
                end
                this.InitialUserData = T(:);
            end
        end
        
        function this = checkInitialErrorValues(this, N)
            % It must be empty or have the same height as X and be all +-1.
            % If empty, create it because there is always this constraint.
            if isempty(this.InitialErrorValues)
                this.InitialErrorValues = NaN(N,1);
            else
                if ~isequal(size(this.InitialErrorValues), [N,1])
                    bayesoptim.err('InitialErrorValuesSize');
                elseif ~all(arrayfun(@bayesoptim.isPlusMinusOneOrNaN, this.InitialErrorValues))
                    bayesoptim.err('InitialErrorValuesType');
                end
                % All checks passed. Convert to double
                this.InitialErrorValues = double(this.InitialErrorValues);
            end
        end
        
        % Hidden public properties
        function this = checkAndFillNumGridDivisions(this)
            if ~(isnumeric(this.NumGridDivisions) && isvector(this.NumGridDivisions))
                bayesoptim.err('NumGridDivisionsType');
            end
            if ~all(isreal(this.NumGridDivisions) & this.NumGridDivisions>1 & ...
                    floor(this.NumGridDivisions)==this.NumGridDivisions)
                bayesoptim.err('NumGridDivisionsType');
            end
            if isscalar(this.NumGridDivisions)
                this.NumGridDivisions = repmat(this.NumGridDivisions, 1, this.VarSpec.NumVars);
            end
            if numel(this.NumGridDivisions) ~= this.VarSpec.NumVars
                bayesoptim.err('NumGridDivisionsSize');
            end
            % Set resolution for categoricals to the number of categories
            isCat = this.VarSpec.isCat;
            if any(isCat)
                NumCats = cellfun(@numel, this.VarSpec.Categories);
                if any(NumCats(isCat) ~= this.NumGridDivisions(isCat))
                    this.NumGridDivisions(isCat) = NumCats(isCat);
                end
            end
        end
        
        function checkNumSeedPoints(this)
%             if ~(bayesoptim.isLowerBoundedIntScalar(this.NumSeedPoints, 1) && ...
%                  isfinite(this.NumSeedPoints))
%                 bayesoptim.err('NumSeedPoints');
%             end
        end
        
        function checkNumInitialTrialsForEachAlgo(this)
            if ~(bayesoptim.isLowerBoundedIntScalar(this.NumInitialTrialsForEachAlgo, 1) && ...
                 isfinite(this.NumInitialTrialsForEachAlgo))
                bayesoptim.err('NumInitialTrialsForEachAlgo');
            end
        end
        
        function checkFitModels(this)
            if ~this.FitModels && ~ismember(this.AcquisitionFunctionName, {'grid','random'})
                bayesoptim.err('FitModels');
            end
        end
        
        function checkAlwaysReportObjectiveErrors(this)
            if ~bayesoptim.isLogicalScalar(this.AlwaysReportObjectiveErrors)
                bayesoptim.err('AlwaysReportObjectiveErrors');
            end
        end
        
        function checkIsClassregRegressionFunction(this)
            if ~bayesoptim.isLogicalScalar(this.IsClassregRegressionFunction)
                bayesoptim.err('IsClassregRegressionFunction');
            end
        end
        
    end
    
    methods(Static)
        function MaxObjectiveEvaluations = DefaultMaxObjectiveEvaluations()
            MaxObjectiveEvaluations = 30;
        end
        
        function MaxTime = DefaultMaxTime()
            MaxTime = Inf;
        end
    end
end

function VarSpec = iCreateVarSpecFromBayesoptVars(BayesoptVars)
% Select vars with Optimize==true
BayesoptVars = BayesoptVars([BayesoptVars.Optimize]);
% Construct VarSpec
VarSpec = bayesoptim.VariableSpecification({BayesoptVars.Name}, {BayesoptVars.Type}, ...
    {BayesoptVars.Transform});
% Set categories and bounds
for v = 1:length(BayesoptVars)
    switch BayesoptVars(v).Type
        case 'categorical'
            strings = BayesoptVars(v).Range;
            catVec = categorical(strings, strings);
            % Store categories, in user-specified order, in categorical
            % type:
            VarSpec.Categories{v} = categorical(categories(catVec), strings);
            VarSpec.LBs(v) = 1;
            VarSpec.UBs(v) = length(VarSpec.Categories{v});
        otherwise % numeric
            if ~all(isfinite(BayesoptVars(v).Range))
                bayesoptim.err('NonfiniteVarRange');
            end
            VarSpec.Categories{v} = [];
            VarSpec.LBs(v) = BayesoptVars(v).Range(1);
            VarSpec.UBs(v) = BayesoptVars(v).Range(2);
    end
end
end

function tf = hasParallel()
tf = logical(exist('parfeval.m','file')) && logical(exist('parfor.m','file'));
end

function Nargout = getObjFcnNargout(ObjFcn, NumCoupledConstraints)
try
    Nargout = nargout(ObjFcn);
catch me
    switch me.identifier
        case 'MATLAB:nargin:isScript'
            % this.ObjectiveFcn is a function defined
            % inside a script.
            return;
        case 'MATLAB:narginout:functionDoesnotExist'
            bayesoptim.err('ObjFcnDoesNotExist', char(ObjFcn));
        otherwise
            rethrow(me);
    end
end
if Nargout == 0
    bayesoptim.err('ObjectiveNargoutZero');
elseif Nargout == 1 && NumCoupledConstraints > 0
    bayesoptim.err('NumCoupledConstraintsMismatch');
end
end

