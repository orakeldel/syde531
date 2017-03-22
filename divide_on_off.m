function [ blower_conf ] = divide_on_off(vect_on_off)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

blower_conf = [];
start_index = 1;
    for i=2:length(vect_on_off)
        if (vect_on_off(i) ~= vect_on_off(i-1)) | (i == length(vect_on_off))
            if i == length(vect_on_off)
                blower_conf = [blower_conf; start_index, i, vect_on_off(i-1)];
            else
                blower_conf = [blower_conf; start_index, i-1, vect_on_off(i-1)];
            end
            start_index = i;
        end
    end
    
end
