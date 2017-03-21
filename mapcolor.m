function y = mapcolor(x)
% map = [0,0,1
%     0.15,0,0.85;
%     0.30,0,0.70;
%     0.45,0,0.55;
%     0.60,0,0.40;
%     0.75,0,0.25;
%     1,0,0;
%     0,0,0];
map = [[28 8 168] ./ 255;
    [47 29 243] ./ 255;
    [57 167 226] ./ 255;    
    [171 183 182] ./ 255;
    [223 143 179] ./ 255;
    [252 53 31] ./ 255;
    [201 25 27] ./ 255;
    [41 10 23] ./ 255;
    ];
    
party = x;
switch party
    case 1
        y = map(1,:);
    case 2
        y = map(2,:);
    case 3
        y = map(3,:);
    case 4
        y = map(4,:);
    case 5
        y = map(5,:);
    case 6
        y = map(6,:);
    case 7
        y = map(7,:);
    otherwise
        y = map(8,:);
end

end