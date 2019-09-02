function [ T ] = StiffnessMatrix(angle)
T = [
    cos(angle)*cos(angle), cos(angle)*sin(angle), -cos(angle)*cos(angle), -cos(angle)*sin(angle);
    cos(angle)*sin(angle), sin(angle)*sin(angle), -cos(angle)*sin(angle), -sin(angle)*sin(angle);
    -cos(angle)*cos(angle), -cos(angle)*sin(angle), cos(angle)*cos(angle), cos(angle)*sin(angle);
    -cos(angle)*sin(angle), -sin(angle)*sin(angle), cos(angle)*sin(angle), sin(angle)*sin(angle)];
end