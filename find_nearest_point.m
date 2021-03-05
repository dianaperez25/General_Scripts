%% this is a script to find point closest to a coordinate
% ------------------------------------------------------------------
function [nearest_point, coordinates, distance] = find_nearest_point(target_coordinates, ordinates)
comp_cords = [];
range = 10;
count = 1;
for c = 1:length(ordinates)
     if ordinates(c,1) < target_coordinates(1)+range && ordinates(c,1) > target_coordinates(1)-range
         if ordinates(c,2) < target_coordinates(2)+range && ordinates(c,2) > target_coordinates(2)-range
             if ordinates(c,3) < target_coordinates(3)+range && ordinates(c,3) > target_coordinates(3)-range
                comp_cords(count,1) = c;
                comp_cords(count,2) = sqrt(((ordinates(c,1)-target_coordinates(1))^2)+((ordinates(c,2)-target_coordinates(2))^2)+((ordinates(c,3)-target_coordinates(3))^2));
                count = count +1;
             end
         end
     end
end

[Y, I] = min(comp_cords);
distance = Y(2);
coordinates = ordinates(Y(1),:);
nearest_point = comp_cords(I(2));

end

    


