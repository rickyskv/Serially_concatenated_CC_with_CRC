% generator matrix of the CRC
function [G_CRC] = generator_CRC(k,m,gen_poly)
    G_CRC=zeros(k,k+m);
    G_CRC(1,:)=[gen_poly,zeros(1,k-1)];
    for i_row = 2:k
        G_CRC(i_row,:)=circshift(G_CRC(1,:),i_row-1);
    end
    for i_row = 2:k
        for j_row = 1:i_row-1
            if G_CRC(j_row,i_row)==1
                G_CRC(j_row,:)=mod(G_CRC(i_row,:)+G_CRC(j_row,:),2);
            end
        end
    end
    %H_CRC=[G_CRC(:,k+1:end);eye(m)].';
    G_CRC = [eye(k),G_CRC(:,k+1:end)];
    
end