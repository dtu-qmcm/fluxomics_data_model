clear functions

samples = table2array(readtable('../../data/poolsize_samples.csv'));

r = reaction({...
...
... % Glycolysis
... % ---------------------------------------------------------------------
      'G6P (abcdef) <-> F6P (abcdef)';...                 % v1
      'F6P (abcdef) -> FBP (abcdef)';...                  % v2
      'FBP (abcdef) <-> DHAP (cba) + GAP (def)';...       % v3
      'DHAP (abc) <-> GAP (abc)';...                      % v4
      'GAP (abc) <-> PG3 (abc)';...                       % v5
      'PG3 (abc) <-> PEP (abc)';...                       % v6
      'PEP (abc) -> Pyr (abc)';...                        % v7
...
... % Pentose phosphate pathway
... % ---------------------------------------------------------------------
      'G6P (abcdef) -> PG6 (abcdef)';...                  % v8
      'PG6 (abcdef) -> Ru5P (bcdef) + CO2 (a)';...        % v9
      'Ru5P (abcde) <-> X5P (abcde)';...                  % v10
      'Ru5P (abcde) <-> R5P (abcde)';...                  % v11
      'X5P (abcde) <-> GAP (cde) + EC2 (ab)';...          % v12
      'F6P (abcdef) <-> E4P (cdef) + EC2 (ab)';...        % v13
      'S7P (abcdefg) <-> R5P (cdefg) + EC2 (ab)';...      % v14
      'F6P (abcdef) <-> GAP (def) + EC3 (abc)';...        % v15
      'S7P (abcdefg) <-> E4P (defg) + EC3 (abc)';...      % v16
...
... % Entner-Doudoroff pathway
... % ---------------------------------------------------------------------
      'PG6 (abcdef) -> KDPG (abcdef)';...                 % v17
      'KDPG (abcdef) -> Pyr (abc) + GAP (def)';...        % v18
...
... % Tricarboxylic acid cycle
... % ---------------------------------------------------------------------
      'Pyr (abc) -> AcCoA (bc) + CO2 (a)';...             % v19
      'OAA (abcd) + AcCoA (ef) -> Cit (dcbfea)';...       % v20
      'Cit (abcdef) <-> ICit (abcdef)';...                % v21
      'ICit (abcdef) <-> AKG (abcde) + CO2 (f)';...       % v22
      'AKG (abcde) -> SucCoA (bcde) + CO2 (a)';...        % v23
      'SucCoA (abcd) <-> Suc (abcd)';...                  % v24
      'Suc (abcd) <-> Fum (abcd)';...                     % v25
      'Fum (abcd) <-> Mal (abcd)';...                     % v26
      'Mal (abcd) <-> OAA (abcd)';...                     % v27
...
... % Amphibolic reactions
... % ---------------------------------------------------------------------
      'Mal (abcd) -> Pyr (abc) + CO2 (d)';...             % v28
      'PEP (abc) + CO2 (d) <-> OAA (abcd)';...            % v29
...
... % Acetic acid formation
... % ---------------------------------------------------------------------
      'AcCoA (ab) <-> Ac (ab)';...                        % v30
...
... % PDO biosynthesis
... % ---------------------------------------------------------------------
      'DHAP (abc) <-> Glyc3P (abc)';...                   % v31
      'Glyc3P (abc) -> Glyc (abc)';...                    % v32
      'Glyc (abc) -> HPA (abc)';...                       % v33
      'HPA (abc) -> PDO (abc)';...                        % v34
...
... % Amino acid biosynthesis
... % ---------------------------------------------------------------------
      'AKG (abcde) -> Glu (abcde)';...                                                                                                                        % v35
      'Glu (abcde) -> Gln (abcde)';...                                                                                                                        % v36
      'Glu (abcde) -> Pro (abcde)';...                                                                                                                        % v37
      'Glu (abcde) + CO2 (f) + Gln (ghijk) + Asp (lmno) + AcCoA (pq) -> Arg (abcdef) + AKG (ghijk) + Fum (lmno) + Ac (pq)';...                                % v38
      'OAA (abcd) + Glu (efghi) -> Asp (abcd) + AKG (efghi)';...                                                                                              % v39
      'Asp (abcd) -> Asn (abcd)';...                                                                                                                          % v40
      'Pyr (abc) + Glu (defgh) -> Ala (abc) + AKG (defgh)';...                                                                                                % v41
      'PG3 (abc) + Glu (defgh) -> Ser (abc) + AKG (defgh)';...                                                                                                % v42
      'Ser (abc) <-> Gly (ab) + MEETHF (c)';...                                                                                                               % v43
      'Gly (ab) <-> CO2 (a) + MEETHF (b)';...                                                                                                                 % v44
      'Thr (abcd) -> Gly (ab) + AcCoA (cd)';...                                                                                                               % v45
      'Ser (abc) + AcCoA (de) -> Cys (abc) + Ac (de)';...                                                                                                     % v46
      'Asp (abcd) + Pyr (efg) + Glu (hijkl) + SucCoA (mnop) -> LL_DAP (abcdgfe) + AKG (hijkl) + Suc (mnop)';...                                               % v47
      'LL_DAP (abcdefg) -> Lys (abcdef) + CO2 (g)';...                                                                                                        % v48
      'Asp (abcd) -> Thr (abcd)';...                                                                                                                          % v49
      'Asp (abcd) + METHF (e) + Cys (fgh) + SucCoA (ijkl) -> Met (abcde) + Pyr (fgh) + Suc (ijkl)';...                                                        % v50
      'Pyr (abc) + Pyr (def) + Glu (ghijk) -> Val (abcef) + CO2 (d) + AKG (ghijk)';...                                                                        % v51
      'AcCoA (ab) + Pyr (cde) + Pyr (fgh) + Glu (ijklm) -> Leu (abdghe) + CO2 (c) + CO2 (f) + AKG (ijklm)';...                                                % v52
      'Thr (abcd) + Pyr (efg) + Glu (hijkl) -> Ile (abfcdg) + CO2 (e) + AKG (hijkl)';...                                                                      % v53
      'PEP (abc) + PEP (def) + E4P (ghij) + Glu (klmno) -> Phe (abcefghij) + CO2 (d) + AKG (klmno)';...                                                       % v54
      'PEP (abc) + PEP (def) + E4P (ghij) + Glu (klmno) -> Tyr (abcefghij) + CO2 (d) + AKG (klmno)';...                                                       % v55
      'Ser (abc) + R5P (defgh) + PEP (ijk) + E4P (lmno) + PEP (pqr) + Gln (stuvw) -> Trp (abcedklmnoj) + CO2 (i) + GAP (fgh) + Pyr (pqr) + Glu (stuvw)';...   % v56
      'R5P (abcde) + FTHF (f) + Gln (ghijk) + Asp (lmno) -> His (edcbaf) + AKG (ghijk) + Fum (lmno)';...                                                      % v57
...
... % One carbon metabolism
... % ---------------------------------------------------------------------
      'MEETHF (a) -> METHF (a)';...                       % v58
      'MEETHF (a) -> FTHF (a)';...                        % v59
...
... % Transport
... % ---------------------------------------------------------------------
      'Gluc.pre (abcdef) -> G6P (abcdef)';...                   % v60
      'Gluc.ext (abcdef) -> G6P (abcdef)';...                   % v61
      'Cit.ext (abcdef) -> Cit (abcdef)';...                    % v62
      'Glyc.ext (abc) + Dummy.ext <-> Glyc (abc) + Dummy';...   % v63
      'PDO (abc) -> PDO.ext (abc)';...                          % v64
      'Ac (ab) -> Ac.ext (ab)';...                              % v65
      'CO2 (a) -> CO2.ext (a)';...                              % v66
...
... % Biomass formation
... % ---------------------------------------------------------------------
      '0.488 Ala + 0.281 Arg + 0.229 Asn + 0.229 Asp + 0.087 Cys + 0.250 Glu + 0.250 Gln + 0.582 Gly + 0.090 His + 0.276 Ile + 0.428 Leu + 0.326 Lys + 0.146 Met + 0.176 Phe + 0.210 Pro + 0.205 Ser + 0.241 Thr + 0.054 Trp + 0.131 Tyr + 0.402 Val + 0.205 G6P + 0.071 F6P + 0.754 R5P + 0.129 GAP + 0.619 PG3 + 0.051 PEP + 0.083 Pyr + 2.510 AcCoA + 0.087 AKG + 0.340 OAA + 0.443 MEETHF -> 39.68 Biomass';...   % v67
...
... % Dummy fluxes
... % ---------------------------------------------------------------------
      'Dummy -> Dummy.ext';...                                  % v68
...
});

% tracers
t = tracer({...
'1-13C_Gluc: Gluc.ext @ 1';...
'U-13C_Gluc: Gluc.ext @ 1 2 3 4 5 6';...
});
t.frac = [0.75,0.25];

% flux measurements
f = data('R61 R68 R64 R66 R67 R62 R65');
f.val = [...
      92.5864,...       % Gluc.ext -> G6P
      1.983114,...      % Glyc.ext -> Glyc 
      129.16588,...     % PDO -> PDO.ext
      177.54056,...     % CO2 -> CO2.ext
      0.90369627,...    % Biomass formation
      0.26614996,...    % Cit.ext -> Cit
      0.31357001...     % Ac -> Ac.ext
];
f.std = f.val/20;
f.val = normrnd(f.val,f.std); % introduce random errors to flux measurements

% MS data
d = msdata({...
      'AKG5: AKG @ 1 2 3 4 5';...
      'Ala2: Ala @ 2 3';...
      'Ala3: Ala @ 1 2 3';...
      'Asp2a: Asp @ 1 2';...
      'Asp2b: Asp @ 1 2';...
      'Asp3: Asp @ 2 3 4';...
      'Asp4: Asp @ 1 2 3 4';...
      'Cit6: Cit @ 1 2 3 4 5 6';...
      'Glu4: Glu @ 2 3 4 5';...
      'Glu5: Glu @ 1 2 3 4 5';...
      'Gly1: Gly @ 2';...
      'Gly2: Gly @ 1 2';...
      'Ile5a: Ile @ 2 3 4 5 6';...
      'Ile5b: Ile @ 2 3 4 5 6';...
      'Leu5: Leu @ 2 3 4 5 6';...
      'Mal4: Mal @ 1 2 3 4';...
      'Met4a: Met @ 2 3 4 5';...
      'Met4b: Met @ 2 3 4 5';...
      'Met5: Met @ 1 2 3 4 5';...
      'Phe8a: Phe @ 2 3 4 5 6 7 8 9';...
      'Phe2: Phe @ 1 2';...
      'Phe8b: Phe @ 2 3 4 5 6 7 8 9';...
      'Phe9: Phe @ 1 2 3 4 5 6 7 8 9';...
      'Pyr3: Pyr @ 1 2 3';...
      'Ser2a: Ser @ 2 3';...
      'Ser2b: Ser @ 1 2';...
      'Ser2c: Ser @ 2 3';...
      'Suc4: Suc @ 1 2 3 4';...
      'Thr3: Thr @ 2 3 4';...
      'Thr4: Thr @ 1 2 3 4';...
      'Tyr2: Tyr @ 1 2';...
      'Val4: Val @ 2 3 4 5';...
      'Val5: Val @ 1 2 3 4 5';...
});
d.idvs = idv; % create devault IDVs, which will be replaced by simulated IDVs later.

% set up model
x = experiment(t);
x.data_flx = f;
x.data_ms = d;
m = model(r);
m.expts = x;

% Take care of symmetrical metabolites
m.mets{'Suc'}.sym = list('rotate180',atommap('1:4 2:3 3:2 4:1'));
m.mets{'Fum'}.sym = list('rotate180',atommap('1:4 2:3 3:2 4:1'));

%% 

% simulate MS measurements
nmts = 8;                               % number of total measurements
samp = 8/60/60;                         % spacing between measurements in hours
% INST
m.options.int_tspan = 0:samp:(samp*nmts);   % time points in hours
% STAT
%m.options.int_tspan = [];   % time points in hours
m.options.sim_tunit = 'h';              % hours are unit of time
m.options.fit_reinit = true;
m.options.sim_ss = false;
m.options.sim_sens = false;

m.options.int_reltol = 10^-6;

m.states{'Glyc.ext'}.bal = false;
m.rates.flx.lb = 1e-7;
m.rates.flx.ub = 1e7;
m.rates.flx.val = [...
      77.5311428582956,...     % G6P -> F6P
      17.1878684051923,...     % F6P -> G6P
      85.2342527495551,...     % F6P -> FBP
      290.460568690959,...     % FBP -> DHAP + GAP
      205.226315941404,...     % DHAP + GAP -> FBP
      274.58349210525,...      % DHAP -> GAP
      316.531346047879,...     % GAP -> DHAP
      22452.7327332828,...     % GAP -> PG3
      22397.1392683668,...     % PG3 -> GAP
      989.267498573664,...     % PG3 -> PEP
      935.246143767407,...     % PEP -> 3PG
      48.3879065275648,...     % PEP -> Pyr
      38.4677367634062,...     % G6P -> PG6
      38.40734676431,...       % PG6 -> Ru5P + CO2
      109.993962875841,...     % Ru5P -> X5P
      85.038822189774,...      % X5P -> Ru5P
      26363.1226030407,...     % Ru5P -> R5P
      26349.6703969624,...     % R5P -> Ru5P
      54.0193540016333,...     % X5P -> GAP + EC2
      29.0642133155664,...     % GAP + EC2 -> X5P
      5.06020365280151,...     % F6P -> E4P + EC2
      17.374656934912,...      % E4P + EC2 -> F6P
      7.00268470770259e-007,...% S7P -> R5P + EC2
      12.6406881042248,...     % R5P + EC2 -> S7P
      2721.0906560371,...      % F6P -> GAP + EC3
      2733.73134344106,...     % GAP + EC3 -> F6P
      17236.9753436607,...     % S7P -> E4P + EC3
      17224.3346562568,...     % E4P + EC3 -> S7P
      0.0603899990961896,...   % PG6 -> KDPG
      0.0603899990961896,...   % KDPG -> Pyr + GAP
      48.6774682891174,...     % Pyr -> AcCoA + CO2
      45.7121925483922,...     % OAA + AcCoA -> Cit
      327.201452172649,...     % Cit -> ICit
      281.223109200753,...     % ICit -> Cit
      459.824613244887,...     % ICit -> AKG + CO2
      413.84627027299,...      % AKG + CO2 -> ICit
      45.004159084501,...      % AKG -> SucCoA + CO2
      46.1578149647674,...     % SucCoA -> Suc
      1.58020021686256,...     % Suc -> SucCoA
      332.33436213708,...      % Suc -> Fum
      287.330203052579,...     % Fum -> Suc
      169779.414115059,...     % Fum -> Mal
      169734.074684896,...     % Mal -> Fum
      249.222338605285,...     % Mal -> OAA
      206.49151079464,...      % OAA -> Mal
      2.60860235198545,...     % Mal -> Pyr + CO2
      11.6730936678694,...     % PEP + CO2 -> OAA
      6.73820210991639,...     % OAA -> PEP + CO2
      1180.31578504477,...     % AcCoA -> Ac
      1180.46671495906,...     % Ac -> AcCoA
      231.285023489945,...     % DHAP -> Glyc3P
      104.102916797761,...     % Glyc3P -> DHAP
      127.182106692184,...     % Glyc3P -> Glyc
      129.16522116879,...      % Glyc -> HPA
      129.16522116879,...      % HPA -> PDO
      5.90235643259762,...     % AKG -> Glu
      0.609994549157668,...    % Glu -> Gln
      0.189776081960163,...    % Glu -> Pro
      0.253938471575266,...    % Glu + CO2 + Gln + Asp + AcCoA -> Arg + AKG + Fum + Ac
      1.64627030655621,...     % Oac + Glu -> Asp + AKG
      0.206946298899416,...    % Asp -> Asn
      0.441003466650284,...    % Pyr + Glu -> Ala + AKG
      1.01272251577938,...     % PG3 + Glu -> Ser + AKG
      1.41211616160504,...     % Ser -> Gly + MEETHF
      0.84401189498947,...     % Gly + MEETHF -> Ser
      0.215792532527811,...    % Gly -> CO2 + MEETHF
      0.170287467472186,...    % CO2 + MEETHF -> Gly
      0.00335165415822219,...  % Thr -> Gly + AcCoA
      0.210561081412943,...    % Ser + AcCoA -> Cys + Ac
      0.294604774852444,...    % Asp + Pyr + Glu + SucCoA -> LL_DAP + AKG + Suc
      0.294604774852444,...    % LL_DAP -> Lys + CO2
      0.470562294031577,...    % Asp -> Thr
      0.131939561743733,...    % Asp + METHF + Cys + SucCoA -> Met + Pyr + Suc
      0.363285642609455,...    % Pyr + Pyr + Glu -> Val + CO2 + AKG
      0.386781728947381,...    % AcCoA + Pyr + Pyr + Glu -> Leu + CO2 + CO2 + AKG
      0.249419993433358,...    % Thr + Pyr + Glu -> Ile + CO2 + AKG
      0.159050430595184,...    % PEP + PEP + E4P + Glu -> Phe + CO2 + AKG
      0.118384127318007,...    % PEP + PEP + E4P + Glu -> Tyr + CO2 + AKG
      0.0487995639326134,...   % Ser + R5P + PEP + E4P + PEP + Gln -> Trp + CO2 + GAP + Pyr + Glu
      0.0813326065543557,...   % R5P + FTHF + Gln + Asp -> His + AKG + Fu
      0.131939561743733,...    % MEETHF -> METHF
      0.0813326065543557,...   % MEETHF -> FTHF
      6.40840700322167,...     % Gluc.pre -> G6P
      92.5878618171061,...     % Gluc.ext -> G6P
      0.266150423504335,...    % Cit.ext -> Cit
      91.9608073960771,...     % Glyc.ext + Dummy.ext -> Glyc + Dummy
      89.9776929194707,...     % Glyc + Dummy -> Glyc.ext + Dummy.ext
      129.16522116879,...      % PDO -> PDO.ext
      0.313569638697864,...    % Ac -> Ac.ext
      177.539702487974,...     % CO2 -> CO2.ext
      0.90369562838173,...     % Biomass formation
      1.98311447660645...      % Dummy -> Dummy.ext
];
m.rates.flx.val = mod2stoich(m)'; % make sure the fluxes are feasible

m.states.lb = 1e-6;
m.states.ub = 1;
m.states.val = [...
      8.41e-5,...       % AKG           1
      1e-5,...          % Ac            2
      0.1,...           % Ac.ext        does not change  3
      1.80e-4,...       % AcCoA         4
      2.10e-3,...       % Ala           5
      1.80e-3,...       % Arg           6
      1.15e-3,...       % Asn           7
      1.38e-3,...       % Asp           8
      0.1,...           % Biomass       does not change  9
      1e-5,...          % CO2           10
      0.1,...           % CO2.ext       does not change 11
      2.68e-3,...       % Cit           12
      0.1,...           % Cit.ext       does not change 13
      1e-4,...          % Cys           14
      1e-4,...          % DHAP          15
      0.1,...           % Dummy         does not change 16
      1,...             % Dummy.ext     does not change 17
      5.87e-5,...       % E4P           18
      1e-4,...          % EC2           19
      1e-4,...          % EC3           20
      3.59e-4,...       % F6P           21
      1.63e-4,...       % FBP           22
      1e-5,...          % FTHF          23
      5.4e-5,...        % Fum           24
      2.08e-3,...       % G6P           25
      1.31e-4,...       % GAP           26
      5.43e-3,...       % Gln           27
      1.38e-2,...       % Glu           28
      0.1,...           % Gluc.ext      does not change 29
      1,...             % Gluc.pre      does not change 30
      1e-4,...          % Gly           31
      1e-4,...          % Glyc          32
      0.1,...           % Glyc.ext      does not change 33
      1e-5,...          % Glyc3P        34
      1e-5,...          % HPA           35
      1.38e-3,...       % His           36
      2.68e-3,...       % ICit          37
      2.18e-4,...       % Ile           38
      1e-5,...          % KDPG          39
      1e-5,...          % LL_DAP        40
      1.34e-3,...       % Leu           41
      8.80e-4,...       % Lys           42
      1e-5,...          % MEETHF        43
      1e-5,...          % METHF         44
      2.79e-4,...       % Mal           45
      1e-4,...          % Met           46
      1e-5,...          % OAA           47
      1e-3,...          % PDO           48
      0.1,...           % PDO.ext       does not change 49
      1.60e-3,...       % PEP           50
      1.28e-3,...       % PG3           51
      4.84e-6,...       % PG6           52
      1.22e-4,...       % Phe           53
      3.02e-4,...       % Pro           54
      1.60e-3,...       % Pyr           55
      2.38e-4,...       % R5P           56
      6.65e-5,...       % Ru5P          57
      1.65e-4,...       % S7P           58
      3.17e-4,...       % Ser           59
      5.47e-5,...       % Suc           60
      1e-5,...          % SucCoA        61
      4.21e-4,...       % Thr           62
      1e-4,...          % Trp           63
      2.60e-4,...       % Tyr           64
      7.65e-4,...       % Val           65
      8.26e-5,...       % X5P           66
];

% order of pools in csv file
indices = [60, 54, 20, 51, 5, 35, 53, 65, 55, 1, 36, 23, 52, 27, 14, 18, ...
    25, 48, 58, 37, 45, 22, 40, 31, 10, 24, 62, 38, 2, 26, 8, 42, 64, ...
    12, 19, 56, 43, 32, 39, 57, 7, 44, 59, 50, 47, 21, 41, 28, 46, 4, ...
    66, 34, 63, 15, 6, 61];
 
% simulate measurements
fprintf("Running pre-evaluation\n");
s = evalc('simulate(m)');
s = evalc('simulate(m)');

%indices = [1,2,4:8,10,12,14:15,18:28,31,32,34:48,50:66];

time = zeros(size(samples, 1),1);
REPS = 5;
temp_states = m.states.val;

results = cell(size(samples,1),1);

fprintf('Progress:\n');
fprintf('    0/% 5d\n', size(samples, 1))

for row = 1:size(samples,1)
    for back = 1:12
        fprintf("\b");
    end
    fprintf('% 5d/% 5d\n', row, size(samples, 1));
    accu = zeros(REPS, 1);
    temp_states(indices) = samples(row,:);
    m.states.val = temp_states;
    for i = 1:REPS
        [s, out] = evalc('simulate(m)');
        a = strsplit(s, "\n");
        b = strsplit(a{4},{':', ' '});
        accu(i) = str2double(b{3});
    end
    results{row} = out.val;
    time(row) = min(accu);
end
disp("======================")
fprintf("runtime %f ms +- %f\n\n", 1000*mean(time), 1000*std(time));

fprintf("%d, ", 1000*time);
fprintf("\n");

save("ecoli_ranges_inca.mat", "results");

