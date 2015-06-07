function [ output_args ] = AISTimeSerie( dataset, training_percentage, ...
    population_percentage, max_iterations, r_threshold )
%AISTimeSerie This is a function for forecasting time series using an
%artificial immune system
%   The function is implemented using the algorithm presented in [2].
%   The proposed AIS (Artificial Immune System) contains immune memory
%   consisting of two populations of Antibodies (ABs). The population of
%   x-antibodies (ABx) detects antigens (AG) representing patterns 
%   x = [x_1,x_2, ...,x_n]^T - AGx, while the populations of y-antibodies
%   (ABy) detects antigens representing patterns 
%   y = [y_1,y_2,...,y_n]^T - AGy. The vectors x and y form the epitopes of
%   AGs and paratopes of ABs. ABx has the cross-reactivity threshold r
%   defining the AB recogition region. The recognition region is
%   represented by the n-dimensional hypersphere of radius r with center at
%   the point x. Similarly ABy has the recognition region of radius s with
%   center at the point y. The cross-reactivity thresholds are adjusted 
%   individually during training. The recognitions regions contain AGs with
%   similar epitopes [2].
%   AG can be bound to many different ABs of the same type (x or y). The
%   strength of binding (affinity) is dependent on the distance between an
%   epitope and a paratope. AB represents a cluster of similar AGs in the
%   pattern space X or Y. The clusters are overlapped and their sizes
%   depend on the similarity between AGs belonging to them, measured in the
%   both patterns space X and Y. The kth ABx can be writtern as a pair
%   {p_k,r_k}, where p_k = x_k and the kth ABy as {q_k,s_k}, where q_k=y_k.
%   After the two populations of immune memory have been created, the
%   empirical conditional probabiliteis P(ABy_k | ABx_j), j,k = 1,2,---,N,
%   that the ith AGy stimulates (is recognized by) the kth ABy, when the
%   corresponding ith AGx stimulates the jth ABx, are determined. These
%   probabilities are calculated for each pair of ABx and ABy on the basis
%   of recoginition fo the training population of AGs [2].
%   In the forecasting phase the new AGx, representing pattern x*, is 
%   presented to the trained immune memory. The forecasted pattern y paired
%   with x* is calculated as the mean of ABy paratopes weighted by the
%   conditional probabilities and affinities.
%
%   Bibliography
%   [1] Dudek, Grzegorz. "Artificial immune system for short-term electric 
%   load forecasting." Artificial Intelligence and Soft Computing–ICAISC 
%   2008. Springer Berlin Heidelberg, 2008. 1007-1017.
%   [2] Dudek, Grzegorz. "Artificial Immune System for Forecasting Time 
%   Series with Multiple Seasonal Cycles." Transactions on Computational 
%   Collective Intelligence XI. Springer Berlin Heidelberg, 2013. 176-197.

switch nargin
    case 1
        training_percentage = 2/3; 
        population_percentage = 0.5;
        max_iterations = 50;
        r_threshold = 10;
    otherwise
        error ('Too few or too many arguments were entered');
end


% Initiliazation of antigens (AG) popoluation C
% Καθορισμός των προτύπων των αντιγόνων ως σύνολο εκπαίδευσης D_T
%while κάποια συνθήκη τερματισμού είναι αναληθείς do
    %for κάθε πρότυπο αντιγόνου z p ∈ D T do
        %Επιλογή ενός υποσυνόλου ΤΛ για έκθεση στο z p ,σαν πληθυσμός S ≤ C;
        %for για κάθε ΤΛ x i ∈ S do
            %Υπολογισμός την ομοιότητα του αντιγόνου μεταξύ z p , x i ;
        %end
        %Επιλογή ενός υποσυνόλου ΤΛ που έχουν την μεγαλύτερη ομοιότητα
        %αντιγόνων σαν πληθυσμός H ≤ S;
        %Προσαρμογή των ΤΛ με κάποια μέθοδο επιλογής, με βάση την
        %υπολογισμένη ομοιότητα και/ή την ομοιότητα του δικτύου
        %των ΤΛ στο H;
        %Ανανέωση του βαθμού ομοιότητας των ΤΛ στο H;
    %end
%end


% Training (imune memory cretion)
% Loading of the training population of antigens (AG).
% An AGx represnets a single x pattern, and AGy represents a single y
% pattern. Both populations of AGx and AGy are divided in the same way  
% into two subsets - training one and test one. The
% first sequences of the time sereis (typically two thirds of the whole
% time series) are included in the taining sets and the latest sequences are
% included in the tests set. Immune memory is trained using the training set
% and after learning the model is test using test set [1,2]
% x is input and y is output

[AG_training,AG_testing,~,~] = SplitTrainingTestingData(dataset,[], ...
                                training_percentage);
                            
[AGx_train,AGx_test,AGy_train,AGy_test] = SplitTrainingTestingData(AG_training, ...
                                AG_testing, population_percentage);                       

% Generation of the antibody populations(AB)
% An initial antibody population is created by copying all the antigens
% from the training set (AB and AG have the same structure).
% This way of initiliazation prevents inserting antibodies in empty regions
% without antigens. Also the paratopes take the form p_k=x_k, q_k=y_k, 
% k=1,2, ...,N. The number of AGs and ABs of both types is the same as the 
% number of learning patterns [1,2]. 

% Πως βρίσκουμε όμως τον αριθμό των patterns ????. Τα θεωρούμε όσες και οι
% μετρήσεις ;;;

ABx_train = AGx_train;
ABx_test = AGx_test;
ABy_train = AGy_etrain;
ABy_test = AGy_test;

paratope_x = [AGx_train, ABx_train];
paratope_y = [AGy_train, ABy_train];

% Calculation of the affinity of antibodies for antigens
% The affinity measure is based on the euclidean distance between x-chains 
% of antigens and antibodies. 
% So d = sqrt(sum_{i=1}^N[x_AB(i) - x_AG(i)]^2), where x_AB and x_AG
% are the x-chains of the antibody and antigen [1]. 

% Calculation of the cross-reactivity thresholds of x-antibodies. The
% recognition region of the kth ABx should be as large as possible and
% cover only the AGx that satisfies two conditions:
% (i)  their epitops x are similar to the paratope p_k, and 
% (ii) the AGy paired with them have epitopes y similar to the kth ABy
%      paratope q_k
% The measure of similarity of the ith AGx to the kth ABx is an affinity and 
% the formulas are in AIS_STLF08.pdf[2]
% Affinity informs about the degree of membership of the ith AGx to the
% cluster represented by the kth ABx.

similiarity_AGx_to_ABx = zeros(size(AGx_train),size(ABx_train));

% Activated antibody detection and evaluation
% If the affinity of the antibody for the antigen is smaller than or equal
% to the cross-reactivity threshold r, it means that the antigen lies in
% the antibody recoginition region (the antibody is activated by the
% antigen). For this antibody the forecast error (MAPE, which is
% traditionally used in short term forecasting) is calulated by:
% delta = 1/24*(sum_{i=1}^24 abs((y_Ab(i) - y_Ag(i))/y_Ag(i))*100%
% where y_Ab and y_Ag are the y chains of the antibody and antigens
% activating this antibody. If several antigens lie in the antibody
% recogintion region, the error is calculated for each of them. The mean
% error is applied to evaluate the antibody and is minimized in the
% following iterations of the algorithm [1].

% Do until the stop criterion is reached
% The algorithm stops if the maximum number of iterations L is reached [1].

	% Clonal selection
    % Each antibody cases sereting as many clones as many antigens are in
    % tis recoginition region. Thus most clones are generated in the dense
    % clusters of antigens [1].    
    
	% Clone hypermutation
    % The main goal of hypermutations is to improve the diversity of the
    % immuyne system in order to effectively recognize new antiges. The
    % hypermutations is realized as foolows. Each clone of the antibody is
    % shifted towards diffwerents antigen lying in the recognition region
    % of this antibody. The bigger the error delta for the given antigen
    % is, the bigger shift toward this antigen is [1]. The formulas are in 
    % AIS_STLF08.pdf. Να ρωτήσω εγώ τι να χρησιμοποιήσω σαν φορμουλες.
	
    % Antibody affinity calculation
    % The affinity measure is based on the euclidean distance between x-chains 
    % of antigens and antibodies. Same as before
	
    % Activated antibody detection and evaluation
    % Same as before
    
	% Selection of the best antibodies
    % For each antigen from the training set, the set of antibodies
    % activated by this antigen is determined. Only one antibody from this
    % set, with the best evalutation delta, is selected to the next
    % population. So the clonal expansion, unnecessary in this model, is
    % halted. The maximum number of antibodies in the next population is
    % equal to the number of antigens, but the real number of antibodies is
    % usually smaller becasuse the same antibody could be selected by
    % several antigens (it depends on the value of the cross-reactivity
    % threshold r). Outlier, i.e. antigen lying away from other antigens,
    % is represented by the seperate antibody [1].


% Calculation of the cross-reactivity thresholds of y-antibodies

% Calculation of the empirical conditional probabilities P(ABy_k | ABx_j)


% Test
% Forecast determination using y-antibodies, probalities P(ABy_k | ABx_j)
% and affinities
% After learning the antibodies rpresent overlapping clusters of similar
% antigens. In the forecast procedure new antigen having only x-chain is
% presented. The \Omega set of antibodies, activated by this antigen, is
% determimed. The y-chains of these antibodies storage average y-chains of
% antigens from the training set with simlar x-chains [1]. The formulas are in 
% AIS_STLF08.pdf

end

