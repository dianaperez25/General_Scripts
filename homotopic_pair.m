function [left_cifti, left_metric, left_coordinates, right_cifti, right_metric, right_coordinates] = homotopic_pair(vertex, ordinates)

load('/projects/p31161/lateralizationVariants/brain_ind.mat');

%% STEP 2: get coordinates for point of interest
test_point = brain_ind(vertex);
test_coords = ordinates(test_point,:);
target_cords = test_coords;
target_cords(1) = target_cords(1)*-1;

%% STEP 3: find distances for vertices in close proximity
comp_cords = [];
range = 10;
count = 1;
for c = 1:length(ordinates)
     if ordinates(c,1) < target_cords(1)+range && ordinates(c,1) > target_cords(1)-range
         if ordinates(c,2) < target_cords(2)+range && ordinates(c,2) > target_cords(2)-range
             if ordinates(c,3) < target_cords(3)+range && ordinates(c,3) > target_cords(3)-range
                comp_cords(count,1) = c;
                %comp_cords(count).distance = sqrt(((right_cords(1)-target_cords(1))^2)+((right_cords(2)-target_cords(2))^2)+((right_cords(3)-target_cords(3))^2));
                comp_cords(count,2) = sqrt(((ordinates(c,1)-target_cords(1))^2)+((ordinates(c,2)-target_cords(2))^2)+((ordinates(c,3)-target_cords(3))^2));
                count = count +1;
             end
         end
     end
end

%% STEP 4: find nearest point

[nearest_point, coordinates, distance] = find_nearest_point(target_cords, ordinates);
pair = find(brain_ind==nearest_point); 

if isempty(pair)
    left_cifti = vertex;
    right_cifti = 0;
    left_metric = test_point;
    right_metric = 0;
    left_coordinates = test_coords;
    right_coordinates = 0;
else
    if nearest_point > 32492 
        left_cifti = vertex;
        right_cifti = pair;
        left_metric = test_point;
        right_metric = nearest_point;
        left_coordinates = test_coords;
        right_coordinates = coordinates;
    else
        left_cifti = pair;
        right_cifti = vertex;
        left_metric = nearest_point;
        right_metric = test_point;
        left_coordinates = coordinates;
        right_coordinates = test_coords;
    end   
end
end

    


