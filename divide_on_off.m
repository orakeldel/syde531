function [ blower_conf ] = divide_on_off(vect_on_off, nc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% nc = number of cycles a day

cins = 24*3600/nc;

blower_conf = [];

last = 0;
for i=1:size(vect_on_off,2)
    blower_conf = [blower_conf; last, (last + vect_on_off(i)), 1];
    next = last + cins;
    blower_conf = [blower_conf; last + vect_on_off(i), next, 0];
    last = next;
end

blower_conf(:,1:2) = blower_conf(:,1:2)/24/3600;


%start_index = 1;
%    for i=2:length(vect_on_off)
%        if (vect_on_off(i) ~= vect_on_off(i-1)) | (i == length(vect_on_off))
%            if i == length(vect_on_off)
%                blower_conf = [blower_conf; start_index, i, vect_on_off(i-1)];
%            else
%                blower_conf = [blower_conf; start_index, i-1, vect_on_off(i-1)];
%            end
%            start_index = i;
%        end
%    end
%    
end
