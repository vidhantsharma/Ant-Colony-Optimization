clc; clear all; close all
%%
NL = 1000; % no. of grids along x
NW = 1000; % no. of grids along y
num_ants = 2;
theta_range = [-90,90]; % in degree
dTheta = 30; % in degree
grid_size = [50,50]; % in km * km % L * W
saturation_limit = [0.1,10]; % for pheromone
alpha = 1;
beta = 5;
v = 250; % m/s
dt = 2;
d = v*dt;
del_d = 50;
MAX_ITR = 500;
R = 1000; % in m
Q = 100;
del_tau0 = 0.1;
% del_taug0 = 100;
del_taug0 = 0;
F = 0.022;

global uniform_search_count
uniform_search_count = 0;
uniform_search_count_th = 1*num_ants;

%% 
% threat Map
global THREAT_MAP
THREAT_MAP.presence = zeros(NW,NL);
THREAT_MAP.threatCentreG = []; % [row,col] % [W,L]
THREAT_MAP.threatCentreP = []; % [x,y]
num_threats = 3;
i = 0;
while i<num_threats
    r1 = randi(NL);
    r2 = randi(NW);
    [x,y] = find_position_from_grid(r1,r2,NL,NW,grid_size);
    dist_between_threats = [];
    for j = 1:size(THREAT_MAP.threatCentreP,1)
        d_threats = norm([x,y] - [THREAT_MAP.threatCentreP(j,1),THREAT_MAP.threatCentreP(j,2)]);
        dist_between_threats = [dist_between_threats;d_threats];
    end
    if any(dist_between_threats<2*R) || x<R || x>grid_size(1)*1000 - R || y<R || y>grid_size(2)*1000 - R
        ;
    else
        THREAT_MAP.threatCentreG = [THREAT_MAP.threatCentreG;r2,r1];
        THREAT_MAP.threatCentreP = [THREAT_MAP.threatCentreP;x,y];
        threat_grids = threat_creator(x,y,NL,NW,R,del_d,grid_size);
        THREAT_MAP.presence(threat_grids(:,1),threat_grids(:,2)) = 1;
        i = i+1;
    end
end
%%
global TARGET_MAP
TARGET_MAP.presence = zeros(NW,NL);
TARGET_MAP.values = zeros(NW,NL);
TARGET_MAP.gridIdx = []; % [row,col]
num_targets = 18;
i = 0;
while i < num_targets
    r1 = randi(NL);
    r2 = randi(NW);
    if (TARGET_MAP.presence(r2,r1) ~= 1 && THREAT_MAP.presence(r2,r1) ~= 1)
        TARGET_MAP.presence(r2,r1) = 1;
        TARGET_MAP.values(r2,r1) = randi([100,200]);
        TARGET_MAP.gridIdx = [TARGET_MAP.gridIdx;r2,r1]; 
        i = i + 1;
    end
end
TARGET_MAP_TEMP = TARGET_MAP;
%%
% SAMSOA
for i=1:num_ants
    ant(i) = create_ant(NL,NW,num_ants,MAX_ITR,i);
end

% filling values for itr 1
% xs = grid_size(1)*1000*rand(num_ants,1); % in m
% ys = grid_size(2)*1000*rand(num_ants,1); % in m
% thetas = theta_range(1) + (theta_range(2)-theta_range(1))*rand(num_ants,1); % in degree 
% % xyts = [xs,ys,thetas];
% xyts = [xs,ys,90*ones(num_ants,1)];

xs = [50;39950];
ys = [50;10050];
thetas = [30;120];
xyts = [xs,ys,thetas];

for i = 1:size(xyts,1)
    ant(i).state(:,:,1) = xyts;
    [grid_L,grid_W] = find_grid(xyts(i,1),xyts(i,2),NL,NW,grid_size);
    ant(i).visited_grids(grid_W,grid_L) = 1;
end

figure(1)
target_plotter(NL,NW,grid_size);
hold on
threat_plotter(R)
hold on

% 
s = "iteration number ";
% main loop
for itr = 1:MAX_ITR
    display(strcat(s,num2str(itr)))

%     plot(ant(1).state(1,1,itr),ant(1).state(1,2,itr),'k.')
%     hold on
%     plot(ant(2).state(2,1,itr),ant(2).state(2,2,itr),'b.')
%     hold on
%     title(strcat(s,num2str(itr)))
%     Frame(itr) = getframe(gcf);
%     drawnow
    

    % state transition
    for j = 1:num_ants
        choosen_candidate(j) = state_transistion(theta_range,dTheta,d,R,ant(j),NL,NW,grid_size,alpha,beta,num_ants,itr);
    end
    
    for j = 1:num_ants
        ant(j) = information_exchange(ant(j),choosen_candidate,num_ants,itr);
    end
    
    % phermone update
    for j = 1:num_ants
        [cand_grid_pheromone,del_tau_inter_ant,del_tau_within_ant] = local_update_pheromone(ant(j),theta_range,dTheta,d,del_d,Q,NL,NW,grid_size,R,del_tau0,num_ants,itr);
        ant(j).pheromoneMap = ant(j).pheromoneMap - del_tau_within_ant;
        ant(j).pheromoneMap = ant(j).pheromoneMap - del_tau_inter_ant;
        % global phermone update
        ant(j).pheromoneMap = ant(j).pheromoneMap + F*del_taug0;
        ant(j).pheromoneMap = pheromone_saturation(ant(j).pheromoneMap,saturation_limit);
    end

visited_grids_global = zeros(NW,NL);
for i = 1:num_ants
    visited_grids_global = visited_grids_global + ant(i).visited_grids;
end
visited_grids_global(find(visited_grids_global>0)) = 1;
for i = 1:num_ants
    ant(i).visited_grids = visited_grids_global;
end
    
end

% SAMSOA plotting

for j = 1:num_ants
    state_plotter(ant(j),grid_size);
    title(strcat(s,num2str(itr)))
    xlabel("in m")
    ylabel("in m")
end

coverage_rate = find_coverage_rate(ant,num_ants,NW,NL); % in percentage
% videoCreator(Frame);
%%
% %random search
% for i=1:num_ants
%     ant(i) = create_ant(NL,NW,num_ants,MAX_ITR,i);
% end
% 
% xs = [50;49950];
% ys = [50;50];
% thetas = [90;90];
% xyts = [xs,ys,thetas];
% 
% for i = 1:size(xyts,1)
%     ant(i).state(:,:,1) = xyts;
%     [grid_L,grid_W] = find_grid(xyts(i,1),xyts(i,2),NL,NW,grid_size);
%     ant(i).visited_grids(grid_W,grid_L) = 1;
% end
% TARGET_MAP = TARGET_MAP_TEMP;
% figure(2)
% target_plotter(NL,NW,grid_size);
% hold on
% threat_plotter(R)
% hold on
%  
% s = "iteration number (random search) ";
% for itr = 1:MAX_ITR
%     display(strcat(s,num2str(itr)))
% 
% %     plot(ant(1).state(1,1,itr),ant(1).state(1,2,itr),'k.')
% %     hold on
% %     plot(ant(2).state(2,1,itr),ant(2).state(2,2,itr),'b.')
% %     hold on
% %     drawnow
%     
% 
%     % state transition
%     for j = 1:num_ants
%         choosen_candidate(j) = random_search_state_transistion(ant(j),d,R,theta_range,dTheta,grid_size,NL,NW,itr);
%     end
%     
%     for j = 1:num_ants
%         ant(j) = information_exchange(ant(j),choosen_candidate,num_ants,itr);
%     end
%     
% end
% %
% % random search plotting
% for j = 1:num_ants
%     state_plotter(ant(j),grid_size);
%     title(strcat(s,num2str(itr)))
%     xlabel("in m")
%     ylabel("in m")
% end
% coverage_rate = find_coverage_rate(ant,num_ants,NW,NL); % in percentage

%%
% %uniform search
% 
% for i=1:num_ants
%     ant(i) = create_ant(NL,NW,num_ants,MAX_ITR,i);
% end
% 
% xs = [50;49950];
% ys = [50;50];
% thetas = [90;90];
% xyts = [xs,ys,thetas];
% 
% for i = 1:size(xyts,1)
%     ant(i).state(:,:,1) = xyts;
%     [grid_L,grid_W] = find_grid(xyts(i,1),xyts(i,2),NL,NW,grid_size);
%     ant(i).visited_grids(grid_W,grid_L) = 1;
% end
% TARGET_MAP = TARGET_MAP_TEMP;
% figure(3)
% target_plotter(NL,NW,grid_size);
% hold on
% 
% s = "iteration number (uniform search) ";
% for itr = 1:MAX_ITR
%     display(strcat(s,num2str(itr)))
% 
% %     plot(ant(1).state(1,1,itr),ant(1).state(1,2,itr),'k.')
% %     hold on
% %     plot(ant(2).state(2,1,itr),ant(2).state(2,2,itr),'b.')
% %     hold on
% %     drawnow
%     
% 
%     % state transition
%     for j = 1:num_ants
%         choosen_candidate(j) = uniform_search_state_transistion(ant(j),d,R,grid_size,NL,NW,uniform_search_count_th,itr);
%     end
%     
%     for j = 1:num_ants
%         ant(j) = information_exchange(ant(j),choosen_candidate,num_ants,itr);
%     end
%     
% end
% 
% %
% % uniform search plotting
% figure(3)
% for j = 1:num_ants
%     state_plotter(ant(j),grid_size);
%     title(strcat(s,num2str(itr)))
%     xlabel("in m")
%     ylabel("in m")
% end
% coverage_rate = find_coverage_rate(ant,num_ants,NW,NL); % in percentage
%% functions
function ant = create_ant(NL,NW,num_ants,num_itr,idx)
    ant.pheromoneMap = ones(NW,NL);
    ant.idx = idx;
    ant.state = zeros(num_ants,3,num_itr); % (x,y,theta)
    ant.visited_grids = zeros(NW,NL);
    ant.attacked_targets = [];
end

function [grid_L,grid_W] = find_grid(x,y,NL,NW,grid_size)
    L = grid_size(1)*1000; % in m
    W = grid_size(2)*1000; % in m
    delta_x = L/NL; % in m
    delta_y = W/NW; % in m

    grid_L = floor(x/delta_x) + 1;
    grid_W = floor(y/delta_y) + 1;
end

function [x,y] = find_position_from_grid(grid_L,grid_W,NL,NW,grid_size)
    L = grid_size(1)*1000; % in m
    W = grid_size(2)*1000; % in m
    delta_x = L/NL; % in m
    delta_y = W/NW; % in m   

    x = (grid_L-1)*delta_x + delta_x/2;
    y = (grid_W-1)*delta_y + delta_y/2;
end

function [candidate_grids,candidate_positions,candidate_steer] = find_candidate_grids(theta_range,dTheta,d,ant,NL,NW,grid_size,num_ants,itr_num)
    global THREAT_MAP
    for i = 1: num_ants
        candidate_grids_all{i} = []; %[x_grid,y_grid]
        candidate_positions_all{i} = []; %[x_pos,y_pos,theta]
        candidate_steer_all{i} = [];
    end
    for i = 1:num_ants
        xcurr = ant.state(i,1,itr_num);
        ycurr = ant.state(i,2,itr_num);
        thetacurr = ant.state(i,3,itr_num);
        for theta = theta_range(1):dTheta:theta_range(2)            
            xnew = xcurr + 2*d*cosd(theta + thetacurr);
            ynew = ycurr + 2*d*sind(theta + thetacurr);
            [cgrid_L,cgrid_W] = find_grid(xnew,ynew,NL,NW,grid_size);
            if xnew>=grid_size(1)*1000 || xnew<=0 || ynew >=grid_size(2)*1000 || ynew<=0 || THREAT_MAP.presence(cgrid_W,cgrid_L)
                    ;
            else
%                 [cgrid_L,cgrid_W] = find_grid(xnew,ynew,NL,NW,grid_size);
                candidate_grids_all{i} = [candidate_grids_all{i};cgrid_L,cgrid_W];
                candidate_positions_all{i} = [candidate_positions_all{i};xnew,ynew,wrapTo360(theta + ant.state(i,3,itr_num))];
                candidate_steer_all{i} = [candidate_steer_all{i};theta];
            end
        end
        if isempty(candidate_grids_all{i})
            [cgrid_L,cgrid_W] = find_grid(xcurr,ycurr,NL,NW,grid_size);
            candidate_grids_all{i} = [candidate_grids_all{i};cgrid_L,cgrid_W];
            possible_turing_angles = [-90,90];
            random_num = randi([1,2]);
            candidate_positions_all{i} = [candidate_positions_all{i};xcurr,ycurr,wrapTo360(ant.state(i,3,itr_num) + possible_turing_angles(random_num))];
            candidate_steer_all{i} = [candidate_steer_all{i};90];
            fprintf("no feasible grids\n")
        end
    end
    common_rows = [];
    for i = 1:num_ants
        if i == ant.idx
            continue;
        else
            [common_row,~] = intersect(candidate_grids_all{ant.idx},candidate_grids_all{i},'rows');
        end
        common_rows = [common_rows;common_row];
    end
    common_rows = unique(common_rows,'rows');

    candidate_grids = candidate_grids_all{ant.idx};
    [candidate_grids,cand_idx] = unique(candidate_grids,'rows');
    candidate_positions = candidate_positions_all{ant.idx};
    candidate_steer = candidate_steer_all{ant.idx};

%     [candidate_grids,cand_idx] = setdiff(candidate_grids,common_rows,'rows');
    candidate_positions = candidate_positions(cand_idx,:);
    candidate_steer = candidate_steer(cand_idx,:);

end

function choosen_candidate = state_transistion(theta_range,dTheta,d,R,ant,NL,NW,grid_size,alpha,beta,num_ants,itr_num)
global TARGET_MAP
[candidate_grids,candidate_positions,candidate_steer] = find_candidate_grids(theta_range,dTheta,d,ant,NL,NW,grid_size,num_ants,itr_num);
obj_vals = [];
Target_vals = [];
for i = 1:size(candidate_grids,1)
    ant_temp = ant;
    curr_grid_L = candidate_grids(i,1);
    curr_grid_x = candidate_positions(i,1);
    curr_grid_W = candidate_grids(i,2);
    curr_grid_y = candidate_positions(i,2);

    d_target_vec = [];
    v_target_vec = [];
    for j = 1:size(TARGET_MAP.gridIdx,1)
        [Target_x,Target_y] = find_position_from_grid(TARGET_MAP.gridIdx(j,2),TARGET_MAP.gridIdx(j,1),NL,NW,grid_size);
        d_target = norm([curr_grid_x,curr_grid_y] - [Target_x,Target_y]);
        v_target = TARGET_MAP.values(TARGET_MAP.gridIdx(j,1),TARGET_MAP.gridIdx(j,2));
        d_target_vec = [d_target_vec,d_target];
        v_target_vec = [v_target_vec,v_target];
    end
    
     
    if (any(d_target_vec<=2*R))
        v_target_vec_temp = v_target_vec(find(d_target_vec<=2*R));
        Target_vals = [Target_vals,max(v_target_vec_temp)];
    else
        Target_vals = [Target_vals,0];
    end
    tau = ant_temp.pheromoneMap(curr_grid_W,curr_grid_L);
    ant_temp.visited_grids(curr_grid_W,curr_grid_L) = 1;
    eta = sum(ant_temp.visited_grids,'all')/(NL*NW);
    obj_val_temp = (tau^alpha) * (eta^beta);
    obj_vals = [obj_vals,obj_val_temp];
end

if (any(Target_vals~=0))
    [val,idx] = max(Target_vals);
    idx1 = find(Target_vals==val);
    candidate_steer_temp = candidate_steer(idx1);
    [~,idx2] = min(abs(candidate_steer_temp));
    val2 = candidate_steer_temp(idx2);
    if val2~=0
        fprintf("target turn");
    end
    idx3 = find(candidate_steer==val2);
    choosen_candidate.mode = "attack";
    choosen_candidate.val = val;
    % also removing the target found from TARGET_MAP
    idx4 = find(v_target_vec == val);
    L = grid_size(1)*1000;
    W = grid_size(2)*1000;
    delta_x = L/NL;
    delta_y = W/NW;
    plot(TARGET_MAP.gridIdx(idx4,2)*delta_x,TARGET_MAP.gridIdx(idx4,1)*delta_y,'g+')
    hold on
    TARGET_MAP.presence(TARGET_MAP.gridIdx(idx4,2),TARGET_MAP.gridIdx(idx4,1)) = 0;
    TARGET_MAP.values(TARGET_MAP.gridIdx(idx4,2),TARGET_MAP.gridIdx(idx4,1)) = 0;
    TARGET_MAP.gridIdx(idx4,:) = [];
else
    rand_num = rand();
    if rand_num<=0.98
        [val,idx] = max(obj_vals);
        idx1 = find(obj_vals==val);
        candidate_steer_temp = candidate_steer(idx1);
        [~,idx2] = min(abs(candidate_steer_temp));
        val2 = candidate_steer_temp(idx2);
         if val2~=0
            fprintf("normal turn");
        end
        idx3 = find(candidate_steer==val2);
        choosen_candidate.mode = "surveillance";
        choosen_candidate.val = val;
    else
        idx3 = randi([1,length(obj_vals)]);
        val = obj_vals(idx3);
        choosen_candidate.mode = "surveillance";
        choosen_candidate.val = val;
    end
end
if length(idx3) > 1
    r = randi(length(idx3));
    idx_val = idx3(r);
else
    idx_val = idx3;
end
choosen_candidate.grid = candidate_grids(idx_val,:);
choosen_candidate.position = candidate_positions(idx_val,:);
end


function ant_updated = information_exchange(ant,choosen_candidate,num_ants,itr_num)
ant_updated = ant;
curr_ant_idx = ant_updated.idx;
curr_ant_visited_grid = choosen_candidate(curr_ant_idx).grid;
curr_ant_visited_position = choosen_candidate(curr_ant_idx).position;
curr_ant_visited_mode = choosen_candidate(curr_ant_idx).mode;
curr_ant_visited_val = choosen_candidate(curr_ant_idx).val;

ant_updated.visited_grids(curr_ant_visited_grid) = 1;

if curr_ant_visited_mode == "attack"
    ant_updated.attacked_targets = [ant_updated.attacked_targets,curr_ant_visited_val];
end

for i = 1:num_ants
    ant_updated.state(i,:,itr_num+1) = choosen_candidate(i).position;
end

end

function [cand_grid_pheromone,del_tau_inter_ant,del_tau_within_ant] = local_update_pheromone(ant,theta_range,dTheta,d,del_d,Q,NL,NW,grid_size,R,del_tau0,num_ants,itr_num)
curr_ant_idx = ant.idx;
x_curr = ant.state(curr_ant_idx,1,itr_num+1);
y_curr = ant.state(curr_ant_idx,2,itr_num+1);
theta_curr = ant.state(curr_ant_idx,3,itr_num+1);
cand_grid_pheromone = [];
cand_pos_pheromone = [];
% d_val = d;
for theta=0:15:360
    for d_val = 0:del_d:R
        x_cand = x_curr + d_val*cosd(theta + theta_curr);
        y_cand = y_curr + d_val*sind(theta + theta_curr);
        if x_cand>=grid_size(1)*1000 || x_cand<=0 || y_cand >=grid_size(2)*1000 || y_cand<=0
            ;
        else
            [cand_grid_L,cand_grid_W] = find_grid(x_cand,y_cand,NL,NW,grid_size);
            cand_grid_pheromone = [cand_grid_pheromone;cand_grid_W,cand_grid_L];
            cand_pos_pheromone = [cand_pos_pheromone;x_cand,y_cand];
        end
    end
end
[cand_grid_pheromone,unique_rows] = unique(cand_grid_pheromone,'rows');
cand_pos_pheromone = cand_pos_pheromone(unique_rows,:);
del_tau_within_ant = zeros(NW,NL);

for i = 1:size(cand_grid_pheromone,1)
%     [cand_x,cand_y] = find_position_from_grid(cand_grid_pheromone(i,2),cand_grid_pheromone(i,1),NL,NW,grid_size);
    cand_x = cand_pos_pheromone(i,1);
    cand_y = cand_pos_pheromone(i,2);
    d_cand_grid_ant = norm([cand_x,cand_y] - [x_curr,y_curr]);
    if d_cand_grid_ant == 0
        del_tau_within_ant(cand_grid_pheromone(i,1),cand_grid_pheromone(i,2)) = 0.9;
    elseif d_cand_grid_ant <=R
%         del_tau_within_ant(cand_grid_pheromone(i,1),cand_grid_pheromone(i,2)) = round((1/(2*1.414*pi))*exp(((Q/d_cand_grid_ant)^2)/2),2);
        del_tau_within_ant(cand_grid_pheromone(i,1),cand_grid_pheromone(i,2)) = round(del_tau0*((R^4 - d_cand_grid_ant^4)/R^4),2);
%         del_tau_within_ant(cand_grid_pheromone(i,1),cand_grid_pheromone(i,2)) = round(del_tau0*((d^4 - d_cand_grid_ant^4)/d^4),3);
    end
end

del_tau_inter_ant = zeros(NW,NL);
for i = 1:num_ants
    if i == curr_ant_idx
        ;
    else
        x_curr = ant.state(i,1,itr_num+1);
        y_curr = ant.state(i,2,itr_num+1);
        theta_curr = ant.state(i,3,itr_num+1);
        del_tau_inter_ant_temp = find_pheromone_decrement(x_curr,y_curr,theta_curr,theta_range,dTheta,d,del_d,Q,NL,NW,grid_size,R,del_tau0);
        del_tau_inter_ant = del_tau_inter_ant + del_tau_inter_ant_temp;
    end
end

end

function del_tau_inter_ant_temp = find_pheromone_decrement(x_curr,y_curr,theta_curr,theta_range,dTheta,d,del_d,Q,NL,NW,grid_size,R,del_tau0)
    cand_grid_pheromone_neigh_ant = [];
    cand_pos_pheromone_neigh_ant = [];
%     d_val = d;
    for theta=0:15:360
        for d_val = 0:del_d:R
            x_cand = x_curr + d_val*cosd(theta + theta_curr);
            y_cand = y_curr + d_val*sind(theta + theta_curr);
            if x_cand>=grid_size(1)*1000 || x_cand<=0 || y_cand >=grid_size(2)*1000 || y_cand<=0
                ;
            else
                [cand_grid_L,cand_grid_W] = find_grid(x_cand,y_cand,NL,NW,grid_size);
                cand_grid_pheromone_neigh_ant = [cand_grid_pheromone_neigh_ant;cand_grid_W,cand_grid_L];
                cand_pos_pheromone_neigh_ant = [cand_pos_pheromone_neigh_ant;x_cand,y_cand];
            end
        end
    end
    [cand_grid_pheromone_neigh_ant,unique_rows] = unique(cand_grid_pheromone_neigh_ant,'rows');
    cand_pos_pheromone_neigh_ant = cand_pos_pheromone_neigh_ant(unique_rows,:);

    del_tau_inter_ant_temp = zeros(NW,NL);

    for i = 1:size(cand_grid_pheromone_neigh_ant,1)
    %     [cand_x,cand_y] = find_position_from_grid(cand_grid_pheromone(i,2),cand_grid_pheromone(i,1),NL,NW,grid_size);
        cand_x = cand_pos_pheromone_neigh_ant(i,1);
        cand_y = cand_pos_pheromone_neigh_ant(i,2);
        d_cand_grid_ant = norm([cand_x,cand_y] - [x_curr,y_curr]);
        if d_cand_grid_ant == 0
            del_tau_inter_ant_temp(cand_grid_pheromone_neigh_ant(i,1),cand_grid_pheromone_neigh_ant(i,2)) = 0.9;
        elseif d_cand_grid_ant <=R
    %         del_tau_within_ant(cand_grid_pheromone(i,1),cand_grid_pheromone(i,2)) = (1/(2*1.414*pi))*exp(((Q/d_cand_grid_ant)^2)/2);
            del_tau_inter_ant_temp(cand_grid_pheromone_neigh_ant(i,1),cand_grid_pheromone_neigh_ant(i,2)) = round(del_tau0*((R^4 - d_cand_grid_ant^4)/R^4),2);
%               del_tau_inter_ant_temp(cand_grid_pheromone_neigh_ant(i,1),cand_grid_pheromone_neigh_ant(i,2)) = round(del_tau0*((d^4 - d_cand_grid_ant^4)/d^4),3);
        end
    end
end

function pheromone_map = pheromone_saturation(pheromone_map,saturation_limit)
    min_sat = saturation_limit(1);
    max_sat = saturation_limit(2);

    for i = 1:size(pheromone_map,1)
        for j = 1:size(pheromone_map,2)
            val = pheromone_map(i,j);
            if val<min_sat
                pheromone_map(i,j) = min_sat;
            elseif val>max_sat
                pheromone_map(i,j) = max_sat;
            else
                ;
            end
        end
    end
end

function choosen_candidate = random_search_state_transistion(ant,d,R,theta_range,dTheta,grid_size,NL,NW,itr_num)
    global TARGET_MAP
    global THREAT_MAP
    curr_ant_idx = ant.idx;
    xcurr = ant.state(curr_ant_idx,1,itr_num);
    ycurr = ant.state(curr_ant_idx,2,itr_num);
    thetacurr = ant.state(curr_ant_idx,3,itr_num);
    candidate_grids = [];
    candidate_positions = [];
    candidate_steer = [];
    for theta=theta_range(1):dTheta:theta_range(2)
        xnew = xcurr + d*cosd(theta + thetacurr);
        ynew = ycurr + d*sind(theta + thetacurr);
        [cgrid_L,cgrid_W] = find_grid(xnew,ynew,NL,NW,grid_size);
        if xnew>=grid_size(1)*1000 || xnew<=0 || ynew >=grid_size(2)*1000 || ynew<=0 || THREAT_MAP.presence(cgrid_W,cgrid_L)
            ;
        else     
            candidate_grids = [candidate_grids;cgrid_L,cgrid_W];
            candidate_positions = [candidate_positions;xnew,ynew,wrapTo360(theta + thetacurr)];
        end
    end
    if (size(candidate_positions,1) == 0)
        [cgrid_L,cgrid_W] = find_grid(xcurr,ycurr,NL,NW,grid_size);
        candidate_grids = [candidate_grids;cgrid_L,cgrid_W];
        candidate_positions = [candidate_positions;xcurr,ycurr,wrapTo360(thetacurr)];      
    end
    [candidate_grids,cand_idx] = unique(candidate_grids,'rows');
    candidate_positions = candidate_positions(cand_idx,:);
    random_idx = randi([1,size(candidate_grids,1)]);
    choosen_candidate.grid = candidate_grids(random_idx,:);
    choosen_candidate.position = candidate_positions(random_idx,:);
    choosen_candidate.mode = 'random search';

    d_target_vec = [];
    v_target_vec = [];
    for j = 1:size(TARGET_MAP.gridIdx,1)
        [Target_x,Target_y] = find_position_from_grid(TARGET_MAP.gridIdx(j,2),TARGET_MAP.gridIdx(j,1),NL,NW,grid_size);
        d_target = norm(choosen_candidate.position(:,1:2) - [Target_x,Target_y]);
        v_target = TARGET_MAP.values(TARGET_MAP.gridIdx(j,1),TARGET_MAP.gridIdx(j,2));
        d_target_vec = [d_target_vec,d_target];
        v_target_vec = [v_target_vec,v_target];
    end
    
     
    if (any(d_target_vec<=2*R))
        v_target_vec_temp = v_target_vec(find(d_target_vec<=2*R));
        target_val = max(v_target_vec_temp);
        [~,idx3] = find(v_target_vec == target_val);
    else
        target_val = 0;
    end
    if (target_val~=0) 
        choosen_candidate.val = target_val;
        % also removing the target found from TARGET_MAP
        L = grid_size(1)*1000;
        W = grid_size(2)*1000;
        delta_x = L/NL;
        delta_y = W/NW;
        plot(TARGET_MAP.gridIdx(idx3,2)*delta_x,TARGET_MAP.gridIdx(idx3,1)*delta_y,'g+')
        hold on
        TARGET_MAP.presence(TARGET_MAP.gridIdx(idx3,2),TARGET_MAP.gridIdx(idx3,1)) = 0;
        TARGET_MAP.values(TARGET_MAP.gridIdx(idx3,2),TARGET_MAP.gridIdx(idx3,1)) = 0;
        TARGET_MAP.gridIdx(idx3,:) = [];
    else
        choosen_candidate.val = nan;
    end
end

function choosen_candidate = uniform_search_state_transistion(ant,d,R,grid_size,NL,NW,uniform_search_count_th,itr_num)
    global uniform_search_count
    global TARGET_MAP
    curr_ant_idx = ant.idx;
    xcurr = ant.state(curr_ant_idx,1,itr_num);
    ycurr = ant.state(curr_ant_idx,2,itr_num);
    thetacurr = ant.state(curr_ant_idx,3,itr_num);
    
    if curr_ant_idx == 1
        xnew = xcurr + d*cosd(thetacurr + 0);
        ynew = ycurr + d*sind(thetacurr + 0);
        thetanew = thetacurr + 0;
        if ynew >=grid_size(2)*1000
            ynew = ycurr;
            thetanew = wrapTo360(thetacurr - 90);
        end
        if ynew <= 0
            ynew = ycurr;
            thetanew = wrapTo360(thetacurr + 90);
        end
        if xnew>xcurr
            uniform_search_count = uniform_search_count + 1;
        end
        if uniform_search_count > uniform_search_count_th
            uniform_search_count = 0;
            if ycurr >= (grid_size(2)*1000)/2
                thetanew = wrapTo360(thetacurr - 90);
            else
                thetanew = wrapTo360(thetacurr + 90);
            end
        end


        d_target_vec = [];
        v_target_vec = [];
        for j = 1:size(TARGET_MAP.gridIdx,1)
            [Target_x,Target_y] = find_position_from_grid(TARGET_MAP.gridIdx(j,2),TARGET_MAP.gridIdx(j,1),NL,NW,grid_size);
            d_target = norm([xnew,ynew] - [Target_x,Target_y]);
            v_target = TARGET_MAP.values(TARGET_MAP.gridIdx(j,1),TARGET_MAP.gridIdx(j,2));
            d_target_vec = [d_target_vec,d_target];
            v_target_vec = [v_target_vec,v_target];
        end
        
        if (any(d_target_vec<=R))
            v_target_vec_temp = v_target_vec(find(d_target_vec<=2*R));
            target_val = max(v_target_vec_temp);
            [~,idx3] = find(v_target_vec == target_val);
        else
            target_val = 0;
        end
        if (target_val~=0) 
            choosen_candidate.val = target_val;
            % also removing the target found from TARGET_MAP
            L = grid_size(1)*1000;
            W = grid_size(2)*1000;
            delta_x = L/NL;
            delta_y = W/NW;
            plot(TARGET_MAP.gridIdx(idx3,2)*delta_x,TARGET_MAP.gridIdx(idx3,1)*delta_y,'g+')
            hold on
            TARGET_MAP.presence(TARGET_MAP.gridIdx(idx3,2),TARGET_MAP.gridIdx(idx3,1)) = 0;
            TARGET_MAP.values(TARGET_MAP.gridIdx(idx3,2),TARGET_MAP.gridIdx(idx3,1)) = 0;
            TARGET_MAP.gridIdx(idx3,:) = [];
        else
            choosen_candidate.val = nan;
        end


    elseif curr_ant_idx == 2
        xnew = xcurr + d*cosd(thetacurr + 0);
        ynew = ycurr + d*sind(thetacurr + 0);
        thetanew = thetacurr + 0;
        if ynew >=grid_size(2)*1000
            ynew = ycurr;
            thetanew = wrapTo360(thetacurr + 90);
        end
        if ynew <= 0
            ynew = ycurr;
            thetanew = wrapTo360(thetacurr - 90);
        end
        if xnew<xcurr
            uniform_search_count = uniform_search_count + 1;
        end
        if uniform_search_count > uniform_search_count_th
            uniform_search_count = 0;
            if ycurr >= (grid_size(2)*1000)/2
                thetanew = wrapTo360(thetacurr + 90);
            else
                thetanew = wrapTo360(thetacurr - 90);
            end
        end

        d_target_vec = [];
        v_target_vec = [];
        for j = 1:size(TARGET_MAP.gridIdx,1)
            [Target_x,Target_y] = find_position_from_grid(TARGET_MAP.gridIdx(j,2),TARGET_MAP.gridIdx(j,1),NL,NW,grid_size);
            d_target = norm([xnew,ynew] - [Target_x,Target_y]);
            v_target = TARGET_MAP.values(TARGET_MAP.gridIdx(j,1),TARGET_MAP.gridIdx(j,2));
            d_target_vec = [d_target_vec,d_target];
            v_target_vec = [v_target_vec,v_target];
        end
        
        if (any(d_target_vec<=R))
            v_target_vec_temp = v_target_vec(find(d_target_vec<=2*R));
            target_val = max(v_target_vec_temp);
            [~,idx3] = find(v_target_vec == target_val);
        else
            target_val = 0;
        end
        if (target_val~=0)
            choosen_candidate.val = target_val;
            % also removing the target found from TARGET_MAP
            L = grid_size(1)*1000;
            W = grid_size(2)*1000;
            delta_x = L/NL;
            delta_y = W/NW;
            plot(TARGET_MAP.gridIdx(idx3,2)*delta_x,TARGET_MAP.gridIdx(idx3,1)*delta_y,'g+')
            hold on
            TARGET_MAP.presence(TARGET_MAP.gridIdx(idx3,2),TARGET_MAP.gridIdx(idx3,1)) = 0;
            TARGET_MAP.values(TARGET_MAP.gridIdx(idx3,2),TARGET_MAP.gridIdx(idx3,1)) = 0;
            TARGET_MAP.gridIdx(idx3,:) = [];
        else
            choosen_candidate.val = nan;
        end

    else
        xnew = xcurr;
        ynew = ycurr;
        thetanew = thetacurr;
    end
    [cgrid_L,cgrid_W] = find_grid(xnew,ynew,NL,NW,grid_size);
    choosen_candidate.grid = [cgrid_L,cgrid_W];
    choosen_candidate.position = [xnew,ynew,thetanew];
    choosen_candidate.mode = 'uniform search';
end

function coverage_rate = find_coverage_rate(ant,num_ants,NW,NL)
visited_grids_all = zeros(NW,NL);    
    for i = 1:num_ants
       visited_grids_all = visited_grids_all + ant(i).visited_grids; 
    end
    visited_grids_all(find(visited_grids_all>0)) = 1;
    coverage_rate = 100*sum(visited_grids_all,'all')/(NL*NW);
end

function state_plotter(ant,grid_size)
x_vec = ant.state(ant.idx,1,:);
x_vec = reshape(x_vec,size(x_vec,3),1);
y_vec = ant.state(ant.idx,2,:);
y_vec = reshape(y_vec,size(y_vec,3),1);

plot(x_vec,y_vec)
xlim([0,grid_size(1)*1000])
ylim([0,grid_size(2)*1000])
hold on
end

function target_plotter(NL,NW,grid_size)
    global TARGET_MAP
    [r,c] = find(TARGET_MAP.presence(:,:)==1);
    L = grid_size(1)*1000; % in m
    W = grid_size(2)*1000; % in m
    delta_x = L/NL; % in m
    delta_y = W/NW; % in m

    for i = 1:length(r)
        plot(c(i)*delta_x,r(i)*delta_y,'r+')
        hold on
    end
end

function threat_grids = threat_creator(x,y,NL,NW,R,del_d,grid_size)
    threat_grids = [];
    for d_val = 0:del_d:R
        for theta = 0:15:360
            xcand = x + d_val*cosd(theta);
            ycand = y + d_val*sind(theta);
            [grid_L,grid_W] = find_grid(xcand,ycand,NL,NW,grid_size);
            threat_grids = [threat_grids;grid_W,grid_L];
        end
    end
    threat_grids = unique(threat_grids,'rows');
end

function threat_plotter(R)
    global THREAT_MAP
    for i = 1:size(THREAT_MAP.threatCentreP,1)
        viscircles([THREAT_MAP.threatCentreP(i,1),THREAT_MAP.threatCentreP(i,2)],R,'color','k');
        hold on
    end
end

function videoCreator(Frame)
    writerObj = VideoWriter('project.avi');

    open(writerObj);

    for i = 1:length(Frame)
        frame = Frame(i);
        writeVideo(writerObj,frame);
    end
end