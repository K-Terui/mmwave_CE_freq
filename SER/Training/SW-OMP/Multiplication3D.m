function MU_matrix = Multiplication3D(mat1,mat2)

%input
%mat1 : (a,c,Z)
%mat2 : (c,b,Z)

%output
%MU_matrix: (a,b,Z)

Mat1 = permute(mat1,[2 1 4 3]);%(c,a,1,Z)
Mat2 = permute(mat2,[1 4 2 3]);%(c,1,b,Z)
M = sum(Mat1 .* Mat2, 1);%(1,a,b,Z)
MU_matrix = permute(M,[2 3 4 1]);


end