% define neighbouring t,f cells (conn 4 or 8)

function [neighs neighs_sb] = get_neighs(conn,t_ind,f_ind,z)
if conn == 4
    % conn = 4
    neighs = [f_ind(z) t_ind(z)-1;f_ind(z) t_ind(z)+1;f_ind(z)-1 t_ind(z);f_ind(z)+1 t_ind(z)];
    neighs_sb = [4 6 8 2];
elseif conn == 8
    % conn = 8
    neighs = [f_ind(z) t_ind(z)-1;f_ind(z) t_ind(z)+1;f_ind(z)-1 t_ind(z);f_ind(z)+1 t_ind(z);...
        f_ind(z)+1 t_ind(z)-1;f_ind(z)+1 t_ind(z)+1;f_ind(z)-1 t_ind(z)-1;f_ind(z)-1 t_ind(z)+1];
    neighs_sb = [4 6 8 2 1 3 7 9];
else
    disp('%%% Conn matrix invalid, valid options 4 and 8 %%%');
end