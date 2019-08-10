fileID = fopen('data.bin');

N = fread(fileID,1,'int');
x = (fread(fileID,N+1,'double'))';
M = fread(fileID,1,'int');
y = (fread(fileID,M+1,'double'))';
Phi = (fread(fileID,[N+1,M+1],'double'))';
F = (fread(fileID,[N+1,M+1],'double'))';


fclose(fileID);
figure
mesh(max(Phi,-1));

%True Surface
    for i=1:N+1
        for j=1:M+1
          Phitrue(i,j) = .1-(x(i)-x(N/2))*(x(i)-x(N/2))-(y(j)-y(M/2))*(y(j)-y(M/2));
        end
    end
    