% Returns the number of pixels in each degree of visual angle. Meaning that
% for some length x in units degree-of-visual-angle, {x * Deg2Pix()} is
% that length in units pixels.
%
% Code taken and modified from params.m.

function d2p = Deg2Pix()
    DISP_NUM =  0;     % The display id of the desired monitor (0 for default).
    DISP_DIAG = 0.39;  % lab CRT 0.546; % for 21.5 in monitor
    SUBJ_DIST = 0.6;   % distance between subject's eye and monitor (meters)

    scrn_res = Screen('Resolution', DISP_NUM);
    W_WIDTH = scrn_res.width/2;
    W_HEIGHT = scrn_res.height/2;
    
    % Below, I'm assuming that the the screen pixels are square (in aggregate)
    screen_diag_pix = sqrt(W_WIDTH^2+W_HEIGHT^2);
    screen_diag_deg = 2*atand(0.5*DISP_DIAG/SUBJ_DIST);
    
    d2p = screen_diag_pix/screen_diag_deg;
end