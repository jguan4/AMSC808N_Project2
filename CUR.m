function [C,U,R] = CUR(A,c,r,k)
C = column_select(A,c,k);
R = column_select(A',r,k);
U = (R\(C\A)')';
R = R';
end