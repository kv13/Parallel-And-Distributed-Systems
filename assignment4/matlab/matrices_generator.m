
% matrices dimensions. Suppose square matrices for convenience.  
% A(NxM), B(NxM)
n_A = 3e6;
m_A = n_A;

n_B = m_A;
m_B = n_A; 

% approximate number of true elements per row
d = 2;

% approximate number of true elements per row for filter
d_f = 4;

% create sparse matrix A
A = sprand(n_A, m_A, d/n_A) > 0;
% create sparse matrix B
B = sprand(n_B, m_B, d/n_B) > 0;
% create sparse matrix F
F = sprand(n_A, m_B, d_f/n_A) > 0;

%********************************************************************%
fprintf('time for bmm\n');
tic;
C = (A*B) > 0;
toc

fprintf('time for bmm filtered\n');
tic;
C_f = (F.*(A*B)) > 0;
toc
%********************************************************************%

% save matrices
% ----A
[rowA,colA] = find(A);
cooA        = sortrows([rowA,colA]);
nnzA        = size(cooA,1);

% save matrix market banner
fid_A       = fopen('../data/matrix_A.mtx','wt');
fprintf(fid_A,'%%%%MatrixMarket matrix coordinate pattern general\n');
fclose(fid_A);

% save elements
csvwrite('../data/matrix_A.mtx', [n_A, m_A, nnzA], 'append', 'on', 'delimiter', ' ', 'roffset', 0);
csvwrite('../data/matrix_A.mtx', cooA, 'append', 'on', 'delimiter', ' ', 'roffset', 0);

% ----B
[rowB,colB] = find(B);
cooB        = sortrows([rowB,colB]);
nnzB        = size(cooB,1);

% save matrix market banner
fid_B       = fopen('../data/matrix_B.mtx','wt');
fprintf(fid_B, '%%%%MatrixMarket matrix coordinate pattern general\n');
fclose(fid_B);

% save elements
csvwrite('../data/matrix_B.mtx', [n_B, m_B, nnzB], 'append', 'on', 'delimiter', ' ', 'roffset', 0);
csvwrite('../data/matrix_B.mtx', cooB, 'append', 'on', 'delimiter', ' ', 'roffset', 0);

% ----F
[rowF,colF] = find(F);
cooF        = sortrows([rowF,colF]);
nnzF        = size(cooF,1);

% save matrix market banner
fid_F       = fopen('../data/matrix_F.mtx','wt');
fprintf(fid_F, '%%%%MatrixMarket matrix coordinate pattern general\n');
fclose(fid_F);

% save elements
csvwrite('../data/matrix_F.mtx', [n_A, m_B, nnzF], 'append', 'on', 'delimiter', ' ', 'roffset', 0);
csvwrite('../data/matrix_F.mtx', cooF, 'append', 'on', 'delimiter', ' ', 'roffset', 0);

% ----C
[rowC,colC] = find(C);
cooC        = sortrows([rowC,colC]);
nnzC        = size(cooC,1);

% save matrix market banner
fid_C       = fopen('../data/matrix_C.mtx','wt');
fprintf(fid_C, '%%%%MatrixMarket matrix coordinate pattern general\n');
fclose(fid_C);

% save elements
csvwrite('../data/matrix_C.mtx', [n_A, m_B, nnzC], 'append', 'on', 'delimiter', ' ', 'roffset', 0);
csvwrite('../data/matrix_C.mtx', cooC, 'append', 'on', 'delimiter', ' ', 'roffset', 0);

% ----C_f
[rowC_f,colC_f] = find(C_f);
cooC_f          = sortrows([rowC_f,colC_f]);
nnzC_f          = size(cooC_f,1);

% save matrix market banner
fidC_f       = fopen('../data/matrix_C_f.mtx','wt');
fprintf(fidC_f, '%%%%MatrixMarket matrix coordinate pattern general\n');
fclose(fidC_f);

% save elements
csvwrite('../data/matrix_C_f.mtx', [n_A, m_B, nnzC_f], 'append', 'on', 'delimiter', ' ', 'roffset', 0);
csvwrite('../data/matrix_C_f.mtx', cooC_f, 'append', 'on', 'delimiter', ' ', 'roffset', 0);