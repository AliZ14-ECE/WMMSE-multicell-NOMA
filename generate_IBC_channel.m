function [H, D, MS, Cell] = generate_IBC_channel(Num_of_user_in_each_cell, inner_radius, Num_of_cell, minR_ratio,seed, d, t ,plot_flag)
% Generate Num_of_users*Num_of_users interference channel
% Num_of_cell cells, each with Num_of_user_in_each_cell user
rng(seed);
% Cell environment channel
Cell.Ncell      = Num_of_cell;   % number of coordinated cells
Cell.Nintra     = Num_of_user_in_each_cell;   % number of cell-intra users
%Cell.Rcell      = inner_radius*2/sqrt(3); % cell radius   %%??? %% By seeing the Layout.m file, it seems like Distance is the half distance between the centers of two adjasent cells.?? but why sqrt(3)???
Cell.Rcell      = inner_radius; % cell radius   %%??? %% By seeing the Layout.m file, it seems like Distance is the half distance between the centers of two adjasent cells.?? but why sqrt(3)???
cells_layout; % Layout of cells

% If consider cellular environment, generate the corresponding channel matrix
MS = usergenerator(Cell, minR_ratio, d, t); %Rcellmin    = minR_ratio*Cell.Rcell;
[HLarge, ~] = channelsample(MS,Cell);
H=(randn(Num_of_user_in_each_cell, Num_of_cell, Num_of_cell)+sqrt(-1)*randn(Num_of_user_in_each_cell, Num_of_cell, Num_of_cell))/sqrt(2); %% The channels between each user in each cell and each basestation.??
Habs = zeros(Num_of_user_in_each_cell, Num_of_cell, Num_of_cell);
I = zeros(Num_of_user_in_each_cell, Num_of_cell);
%for Base=1:BaseNum
%    for base=1:Num_of_cell
        for c=1:Num_of_cell
            for User=1:Num_of_user_in_each_cell
                H(User, c, c)=H(User, c, c)*sqrt(HLarge(User, c, c));
                %between user in Cell MS and base station base
            end
            
            %sorting the users in descending order of their channels
            [Habs(:, c, c), I(:,c)] = sort(abs(H(:, c, c)), 'descend');
            
            H(:, c, c) = H(I(:,c), c, c);

            MS.Position{c}(:, :) = MS.Position{c}(I(:,c), :);


            %plotiing users (after sorting)
            if plot_flag ==1
            for User = 1:Num_of_user_in_each_cell
              plot(MS.Position{c}(User, 1),MS.Position{c}(User, 2), 'o','MarkerSize',3);
              text(MS.Position{c}(User, 1),MS.Position{c}(User, 2),num2str(User),'FontSize',10);
            end
            end
        end
%    end
%end

%regenerating the channel samples after sorting
[HLarge, D] = channelsample(MS,Cell);
%generating the channels:
for base=1:Num_of_cell
  for c=1:Num_of_cell
    for User=1:Num_of_user_in_each_cell
      if base ~= c
        H(User, c, base)=H(User, c, base)*sqrt(HLarge(User, c, base));
      end
    end
  end
end
      
end


function [Hlarge, D] = channelsample(MS,Cell)
Ncell	= Cell.Ncell;      % # of cells
Nintra  = Cell.Nintra;     % # of MSs in each cell
%Nbs     = Cell.NintraBase; % # of BSs in each cell
BS_position = Cell.Position;
% Channel between BSs and Cell-intra MSs
Hlarge	  = zeros( Nintra, Ncell, Ncell );
D = zeros(Nintra, Ncell, Ncell);
% large-scale fading
%% this computes the distances (then the channel) between each user and each BS??
for base = 1 : Ncell %% this combined with the third loop: for each BS in all cells??
    for CellMS = 1 : Ncell %% this combined with the fourth loop: for each user in all the cells??
            for User = 1 : Nintra
                    d = norm(MS.Position{CellMS}(User,:)-BS_position(base,:));
                    D(User, CellMS, base) = d;
%                     PL= (200/d)^(4); %% ??? 
%                     PL= 10^(randn*8/10)*(200/d)^(3);
                    PL=128.1+37.6*log10(d/1000)+8*randn;
                    PL=10^(-PL/10);
                    Hlarge(User,CellMS,base)=PL;
            end
    end
end
end


function MS = usergenerator(Cell, minR_ratio, d, t)
Ncell       = Cell.Ncell;         % number of cells
Nintra      = Cell.Nintra;        % number of cell-intra users
MS.Position = [];
Nms         = Nintra;
% Rcellmin    = minR_ratio*Cell.Rcell;

% User Deployment
MS.Position	= cell(Ncell,1);

if Nms>=1
    for n = 1 : Ncell % generate users for each cell     
        x = d(:,n).*cos(t(:,n));
        y = d(:,n).*sin(t(:,n));
        MS.Position{n}(:,1)	= x+Cell.Position(n,1);  %%  num_users_in_each_cell × 2   matrex. the raw i is the x-y coordinates of user i.
        MS.Position{n}(:,2) = y+Cell.Position(n,2);
    end
end

end
