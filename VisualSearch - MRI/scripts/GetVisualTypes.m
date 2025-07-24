% Uses the human-legible visual type description to form a new structure
% that can be more directly used with the PsychToolbox functions.
function visualTypes = GetVisualTypes(VisualTypesLegible, window)

global luminance_LUT
global color_LUT

% Using a cellarray, define the features (color, motion, etc.) of all the
% possible visual types. VisualTypes{1} will always be the TARGET's
% features.
% 1     2          3          4         5     6            7
% Name  RGB-Color  Radius-px  Direct.   Freq  Phase-(rad)  Ampl-px

deg2pix = Deg2Pix();  % Deg-to-pixel conversion for degrees of visual angle.


visualTypes = cell(size(VisualTypesLegible));
for i = 1:size(visualTypes,2)
    % Preserve (or scale) Name, Size, Freq, Phase, Ampl... respectively.
    visualTypes{i}{1} = VisualTypesLegible{i}{1};
    visualTypes{i}{3} = VisualTypesLegible{i}{3} * deg2pix;
    visualTypes{i}{5} = VisualTypesLegible{i}{5};
    visualTypes{i}{6} = VisualTypesLegible{i}{6} * 2 * pi;
    visualTypes{i}{7} = VisualTypesLegible{i}{7} * deg2pix;
    
    % Convert color names to RGB values.
    switch VisualTypesLegible{i}{2}
        case 'red'
            visualTypes{i}{2} = color_LUT(5001,:);
        case 'green'
            visualTypes{i}{2} = color_LUT(1,:);
        case 'black'
            visualTypes{i}{2} = luminance_LUT(1,:);
        case 'white'
            visualTypes{i}{2} = luminance_LUT(5001,:);
        otherwise
            disp('Error: Unknown color found in visualTypesLegible.');
    end
    
    % Convert motion direction to xy vectors.
    switch VisualTypesLegible{i}{4}
        case 'none'
            visualTypes{i}{4} = [0,0];
        case 'horz'
            visualTypes{i}{4} = [1,0];
        case 'vert'
            visualTypes{i}{4} = [0,1];
        otherwise
            disp('Error: Unknown motion found in visualTypesLegible.');
    end
    
     % Load a matrix representing the requested shape.
    switch VisualTypesLegible{i}{8}
        case 'circle'
            load('media/circle.mat', 'img');
        case 'square'
            load('media/square.mat', 'img');
        case 'triangle'
            load('media/triangle.mat', 'img');
        otherwise
            error('Error: Unknown shape found in visualTypesLegible.');
    end
    img(:,:,1) = img(:,:,1) * visualTypes{i}{2}(1); % Apply the desired colors
    img(:,:,2) = img(:,:,2) * visualTypes{i}{2}(2);
    img(:,:,3) = img(:,:,3) * visualTypes{i}{2}(3);

    %Create texture, store pointer to it.
    visualTypes{i}{8} = Screen('MakeTexture', window, img);
    
end           
