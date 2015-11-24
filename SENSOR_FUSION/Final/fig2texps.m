%% |fig2texPS| - Export Matlab Figure to LaTeX
%  Last Modification
%  Date:    16.11.2011
%  Author:  Adin Ramirez
%  Version: 0.5.5
%
% Original Author:  Peter Scholz (03.11.2010)
% Email:   contact@peter-scholz.net
%
% This MATLAB function converts a figure to a LaTeX file with PSTricks
% compatible macros. The packages |pstricks|, |pst-node|, |pst-plot|,
% |pst-circ| and |moredefs| are required. A minimal example latex code is
% generated at each time this script is run. The tex file can also be a
% stand-alone file so that the paper size matches the figure outline.

%% Motivation
%
% The conversion tool has been written for the following reasons (compared
% to eps with psfrag):
%
% * With |eps with psfrag| it is cumbersome to replace all tick labels
% * When replacing text in a legend with |eps with psfrag|, the size of the
% legend's box is NOT
% adapted -> this is different in this tool because PSTricks macros are
% used and the box size is adapted to its content
% * When resizing an |eps|-Figure, the widths of all lines are scaled as
% well -> With this tool, all linewidths are equal independent on the size of
% the figure. This guarantees a neat appearance in documents with many
% different figures
% * |eps|-Figures, created by Matlab have a surrounding box that may be
% unwanted sometimes, especially in two-column documents. -> Here, the box
% size can be chosen to be either minimal or fixed
% * The output of this tool is a tex file that is easy to read. For
% instance, each grid line is one line in the TeX document. This offers
% the possibility to finetune the figure with any of the powerful PSTricks
% macros

%% Properties
% The aim of this function is to create a tex-code picture of a Matlab
% figure in PSTricks syntax that can still be adapted after creation.
%
% *Supported Features:*
%
% * 2D plots with arbitrary number of lines and patches
% * Linear and logarithmic axes and second y-axis
% * Grids, subgrids, ticks and subticks
% * Arbitrary colors, linestyles and markers
% * Labels, title, legend and colorbar
% * Annotations (TextArrows and TextBoxes) are supported
% * GUI, where the fontsize, plotsize and the export file location can be
% set
% * Lines with many datapoints can be compressed in order to reduce file
% size and compilation time (copied from |Matfig2PGF| by Paul Wagenaars in
% 2008)
% * Convergence Triangles can be shown for logarithmic plots
% * Support of additional code in PSTricks syntax
%
% *Changes to Version 0.52 from 10.06.2010:*
%
% * small bugfixes
% * dataplot -> psline replaced, option linejoin=1 introduced
% * marker more bold
% * fixed a bug for vertical lines that produced NaN values when cut
% * enhancement for patches which are cut at the axes limits
% * Colorbar enhancement, if the tight inset of the right side is larger,
% this value is used instead
% * fixed a bug for axes without tick labels
%
% *NOT Supported in this version:*
%
% * 3D Plots or contour plots
% * Subfigures (if a figure with multiple axes is chosen, a dialog box
% appears that allows to activate a certain axes)
% * Other annotations than TextArrows, and TextBoxes
% * Reverse Axes
%
% The matlab code is based on the
% |<http://www.mathworks.de/matlabcentral/fileexchange/7884 fig2tex>|
% file by Ercan Solak in 2005.
% The key modifications are that the ratios can be set, implementation of
% grids and logarithmic scale. The algorithm to reduce the number of large
% data sets is copied from
% |<http://www.mathworks.com/matlabcentral/fileexchange/12962 Matfig2PGF>|
% by Paul Wagenaars in 2008.
% Bertrand Scherrer added the support of textboxes on 09.06.2010

%% Bugs and ToDo's
%
% * 26.02.09: The annotations (TextArrows) are not shown in the visualized
% figure, this has to be implemented
% * 20.05.10: If patches are cut, the algorithm works in principle but the
% patches are divided into multiple patches according to the algorithm used
% for lines. This is wrong, instead a merging must be enforced in order to
% draw a straight line at the border...
% * 29.06.10: Bug if a legend is added for patches

%% Syntax
%  fig2texPS()
%  fig2texPS(fignum[,'Option',value,...])
%  fig2texPS(fignum,fname[,'Option',value,...])
%  fig2texPS(fignum,globals[,'Option',value,...])

%% Description
% *|fig2texPS()|* exports either the activated figure (|gcf|) or a test
% function, if no figure is present.
%
% *|fig2texPS(fignum)|* exports the figure with the number |fignum|.
%
% *|fig2texPS(fignum,fname)|* exports the figure with the number |fignum|
% and sets its name to |fname|. The name can be changed in the
% Figure Export Properties later on.
%
% *|fig2texPS(fignum,globals)|* exports the figure with the number
% |fignum|. The struct |globals| can be used in order to set different
% settings from outside this function. See the fieldnames for further
% explanation.

%% Figure Export Properties
% <<GUI_settings.png>>
%
% The _Figure Export Properties_ GUI is opened together with a replot of the
% figure for finetuning the export behavior. Each time a paramter is
% changed, the replot is updated automatically. The GUI can be used to set
% the following paramters:
%
% * *Width / cm* The width of either the figure or the axes can be set in
% centimeters. This parameter can be changed in the tex-file as well. This
% feature allows all figures in a
% document to have the same size and style.
% * *Figure or Axes* This checkbox allows to set if the chosen width
% belongs to the axes or the total figure, in case of Figure, the size of
% the axes is chosen automatically, in case of Axes the space for all
% 'outer' text is added.
% * *Height / Width* The height of the axes can be set by the height factor
% relatively to the axes width. This allows for using the same aspect
% ratios of all figures in a latex document.
% * *Axis equal* Alternatively, the height can be set in such a way that
% the data of the x-axis and the y-axis are in the same aspect ratio.
% * *Square* In Addition, the height can be set in such a way being equal
% to the width.
% * *Auto Slider* This slider allows to adjust the surrounding space by
% increasing or decreasing the size of the figure's text. Note, that the
% fontsize has NO influence on the fontsize in the LaTeX document.
% * *Center axes horizontally* Here it can be decided whether the axes should be
% centered horizontally or not. If set, extra space is added either left or
% right of the axes in order to keep the axes centered.
% * *Manual Space Setting* The four fields *Left*, *Bottom*, *Right* and
% *Top* allow to set the surrounding space manually. This is important,
% when using subfigures in the TeX document. Hereby, it is guaranteed that
% all figures are aligned correctly.
% * *Draw figure outline* If activated, a box is plotted around the latex figure
% in order to visualize the real figure size. If deactivated, no box is
% drawn. However, the box can still be activated in the tex-file by
% uncommenting the coresponding line.
% * *Stand-alone file* If activated, a stand-alone tex-file is created that
% can be compiled with dvi -> ps -> pdf. The borders are matched with the
% size of the figure. Using this method, a pdf file can be generated that
% can be used later on for pdftex documents.
% * *Reduce data points* The number of data points can be reduced for lines
% without any markers if this radiobutton is enabled. The reduction factor
% can be set. If the factor gets smaller, the reduction becomes less. A
% zero factor equals no reduction.
% * *Legend style* The style of the legend can be set. Currently, a
% rectangular box, a box with rounded corners, a shadowbox and no box are
% supported.
% * *Export to file* The filename can be set or changed.


%% Example file
% The following example picture shows the functionality of the script:
%
% <<fig2texPS_example.png>>


function fig2texps(varargin)

EPS=1e-10;
s = warning('query', 'all'); % save warnings
warning('off','all'); % disable all warnings
% clc % clear screen
%% INITIALIZATION

% create new struct 'globals' that contains all plot settings with the
% following fields:
globals = struct ( ...
    'ShowGUI',              1, ...      % determine, whether the GUI is displayed or not
    'BoundingBox',          0, ...      % draw a figure outline box for test aspects
    'plotwidth',            10, ...     % width in cemtimeters (Originally: 10)
    'widthType',            'axes', ... % alternatives: 'figure', 'axes'. When specifying figure,
    ...                                 % a fixed surrounding (tightInset_cm) is not supported
    'labelsep',             10, ...     % separation of labels
    'labelsepMin',          8, ...      % minimum Label separation
    'labelsepMax',          20, ...     % maximum Label separation
    'AxesLineWidth',        0.5, ...    % in pt
    'AxesBoxEnd',           0, ...      % if 1, the axes box will be drawn at the end and may overlay curves,
    ...                                 % if 0, the axes box will be drawn before the content
    'TextArrowLineWidth',   0.5, ...    % in pt
    'TextArrowTextColor',   0, ...      % if 1, the color of the text will be same as the arrow,
    ...                                 % if 0, the color of the text will be black
    'TextBoxLineWidth',     0.5, ...    % in pt
    'TextBoxTextColor',     0, ...      % if 1, the color of the text will be same as the box,
    ...                                 % if 0, the color of the text will be black
    'GridGreyFac',          0.8, ...    % grey factor of grid line color (0 ... black, 1 ... white)
    'GridLineWidth',        0.4, ...    % in pt
    'GridLineDotSep',       0.4, ...    % in pt
    'MinorGridLineWidth',   0.4, ...    % in pt
    'MinorGridLineDotSep',  0.8, ...    % in pt
    'TickMajorLengthRel',   0.012, ...  % relative to x-length
    'TickMinorLengthRel',   0.008, ...
    'MarkerSize',           4, ...      % in pt
    'LineWidth',            0.7, ...    % in pt
    'EdgeLineWidth',        0.4, ...    % in pt (linewidth of the edges of the patches)
    'FontSizeXYlabel',      '\small', ... % size of the xlabel, ylabel and title, may be any
    ...                                 % LaTeX compatibel size like \tiny \small etc.
    'FontSizeTickLabels',   '\footnotesize', ... % size of all labels
    'FontSizeLegend',       '\small', ...
    'FontSizeAnnotations',  '\small', ...
    'FloatingPointFormat',  '%6.6f', ...
    'PointReduceFactor',    0.0001, ... % factor that determines the reducing of the plot points
    'idx_LegendStyle',      3, ...      % 1:    'Rectangular Box'
    ...                                 % 2:    'Rounded corners'
    ...                                 % 3:    'Shadowbox'
    ...                                 % 4:    'No Box'
    'LegendOffsetFac',      2, ...      % Factor relative to the TickMajorLength for the offset of the Legend Box
    'heightFactor',         [], ...     % factor between axes width and height. If
    ...                                 % empty [], ratio is set automatically,
    ...                                 % e.g. 0.7 means that the height has 70 % of the width
    'centerAxes',           0, ...      % if 1, extra space is addd to the right side of the axes in order to center the figure
    'tightInset_cm',        [], ...     % The extra space between the axes and the figure borders can be set in cm
    ...                                 % [left bottom right top] in cm or leave empty []
    ...                                 % the width must be the width of the axes
    'ColorBarWidth_cm',     [0.6, 0.3], ...     % Either empty (auto size), scalar (width) or 1x2 vector (width, offset)
    'pathname',             pwd, ...
    'filename',             'fig1.tex', ... % test tex-file name
    'StandAloneFile',       0 ...      % if set to 1, a standalone file is generated that can be used to create a pdf file, it uses the standalone package for faster builds.
    );

% Optional Convergence Triangle Settings:
globals.TriangleOrder = []; % vector with orders that are shown, if empty, no convergence triangle
globals.TriangleColor = [0.4 0.4 0.4]; % rgb vector
globals.TriangleLxRel = 0.3;
globals.TriangleStr   = {};

globals.AddCode = {}; % additional code before '\end{pspicture}'


%% EVALUATE INPUT ARGUMENTS
numInputs=size(varargin,2); % number of input elements
if numInputs==0,            % test case without any inputs
    globals.StandAloneFile = 1; % do it in standalone for zero inputs
    % test if there is at least one figure open
    if isempty(get(0,'CurrentFigure')), % case where no figure is open, create test figure
        fignum=figure; % create new figure and set fignum to this figure
        time=linspace(0,40,500); % create time vector
        Z1=exp(-0.1*time).*sin(time);
        Z2=exp(-0.1*time);
        Z3=-exp(-0.1*time);
        plot(time,Z1,'-b','Linewidth',2);
        hold on; grid on;
        plot(time,Z2,'--r','Linewidth',2);
        plot(time,Z3,'--r','Linewidth',2);
        set(gca,'YLim',[-1,1]);
        xlabel('Time $t$','interpreter','none')
        ylabel('Magnitude','interpreter','none')
        title('This is a test figure created by fig2texPS','interpreter','none')
        legend({'$e^{-0.1 t}\cdot\sin(t)$','$\pm\,e^{-0.1 t}$'},'Location','EastOutside','interpreter','none')
        annotation(gcf,'textarrow',[0.4625 0.367857142857143],...
            [0.706142857142857 0.642857142857143],'TextEdgeColor','none',...
            'String',{'Envelope'});
    else
        fignum=figure(gcf); % choose current figure
    end
    
elseif numInputs==1,    % figure number is given
    fignum=varargin{1};   % copy first input argument to fignum
    
elseif numInputs >= 2,      % figure number and filename are given
    fignum = varargin{1};   % copy first input argument to fignum
    
    if isstruct (varargin{2}) % option |fig2texPS(fignum,globals)|
        glob_cp = varargin{2};
        names = fieldnames (glob_cp);
        for nm = 1 : length (names),
%              globals = setfield (globals, names{ii}, getfield (glob_cp, names{ii}));
            globals.(names{nm}) = glob_cp.(names{nm});
        end
        if numInputs > 2, parse(varargin{3:end}); end; % parse the rest of the inputs, option |fig2texps(fignum,globals,'Option',val,...)|           
    elseif isfield(globals,varargin{2}) % option |fig2texps(fignum,'Option',val,...)|
        parse(varargin{2:end}); % parse the rest of the inputs, 
    else % option |fig2texps(fignum,filename,'Option',val,...)|
        fname  = varargin{2};   % copy second input argument to fname
        [pathstr, name, ~] = fileparts(fname);
        if ~isempty(pathstr),
            globals.pathname = pathstr;
        end
        globals.filename = strcat(strrep(name,'.tex',''),'.tex'); % append .tex if not already there
        
        if numInputs > 2, parse(varargin{3:end}); end; % parse the rest of the inputs, option |fig2texps(fignum,globals,'Option',val,...)|
    end        
else
    h=warndlg('Wrong input setup!','Warning!','modal');
    uiwait(h);
    warning(s); % turn on specified warnings again
    return
end

% test if this figure exists
fig_nums=get(0,'Children');
if isempty(find(fig_nums==fignum, 1))
    h=warndlg('The chosen figure does not exist!','Warning!','modal');
    uiwait(h);
    warning(s); % turn on specified warnings again
    return
end

figure(fignum); % activate the right figure
fh = gcf; % figure handle
set(fh,'WindowStyle','normal')

% get the number and order of axes (subfigures)
childs=get(fh,'children'); % all axes and legends, but wrong order and doubled
subAxes=[]; % init
for ch=1:length(childs), % scan all childs and search for legends
    if strcmp(get(childs(ch),'type'),'legend') || strcmp(get(childs(ch),'type'),'axes')
        legh = legend(childs(ch));
        if isempty(legh) || childs(ch)~=legh, % new axes found
            subAxes=[subAxes, childs(ch)];
        end
    end
end
subAxes=fliplr(subAxes); % extract subaxes by flipping them
numAxes=length(subAxes);
if numAxes==0, % no axis, cancel
    warndlg('No Axes found!')
    warning(s); % turn on specified warnings again
    return
end
ahSecond = []; % init
ColBar   = [];
if numAxes==1, % only one axes, take it
    ah=subAxes(1);
else % multiple axes found
    % check if two axes belong together (second y-axis)
    pos_last=[0 0 0 0]; % init
    numSubplots=0; % init
    for axs=1:numAxes,
        
        % test for colorbars
        if strcmp(get(subAxes(axs),'Tag'), 'Colorbar')            
            
            % save colorbar for later usage and do not treat as a normal
            % axes
            ColBar = get(subAxes(axs));          
            ColBar.cmap = colormap; % get the corresponding colormap
            
        else
            
            set(subAxes(axs),'Units','normalized');
            pos_act=get(subAxes(axs),'Position');
            if norm(pos_act-pos_last)<EPS, % second axes found
                disp('Found two axes that belong together')
                AxesIndxTot(axs)=axs-1;
            elseif strcmp(get(subAxes(axs),'YAxisLocation'),'right') % axes with right y-axis location
                disp('Found axes with right y-axis that has nonmatching size. Converting size ...')
                %             set(subAxes(axs),'Position',pos_last);
                AxesIndxTot(axs)=axs-1;
            else
                numSubplots=numSubplots+1;
                % extract y-label as titles for the subfigures
                Strs{numSubplots}=get(get(subAxes(axs),'YLabel'),'String');
                if isempty(Strs{numSubplots}) % no ylabel
                    Strs{numSubplots}=['Axis ',num2str(numSubplots)]; % generate one
                end
                AxesIndxTot(axs)=axs;
                AxesIndx(numSubplots)=axs;
            end
            pos_last=pos_act; % update
            
        end
    end
    if numSubplots>1, % choose
        [selectionNumber,ok] = listdlg('PromptString','Choose Axis:','ListSize',[160 160],...
            'SelectionMode','single',...
            'ListString',Strs);
    else
        selectionNumber=1; ok=1;
    end
    if not(ok), 
        warning(s); % turn on specified warnings again
        return,
    else
        ah=subAxes(AxesIndx(selectionNumber)); % take chosen axis
        idx=find(AxesIndxTot==AxesIndx(selectionNumber));
        if length(idx)==2; % copy second axes
            ahSecond=subAxes(idx(end));
        end
    end
end

% create empty string cell array in order to rearrange lines later on
StrFile={};

%% EXTRACT ALL IMPORTANT FIGURE DATA

% Test if one or more axes have reverse style -> warning
if strcmp(get(ah,'XDir'),'reverse') || strcmp(get(ah,'YDir'),'reverse') ...
        || strcmp(get(ahSecond,'YDir'),'reverse')
    h=warndlg('At least one axis is reversed! This is currently not supported. If you are using the spy-command, simply flipud(A) and change the labels manually.','Warning!','modal');
    uiwait(h);
end

% Test if Scales are log
if strcmp(get(ah,'XScale'),'log'), fig.xLog=1; else fig.xLog=0; end
if strcmp(get(ah,'YScale'),'log'), fig.yLog=1; else fig.yLog=0; end

% Test if axes are visible
if strcmp(get(ah,'Visible'),'off')
    fig.AxesVisible = 'off';
else
    fig.AxesVisible = 'on';
end

% Try to extract bases by comparing Ticks with TickLabels
XBase   = [];
xTicks  = get (ah, 'XTick');
xLabels = get (ah, 'XTickLabel'); % array of chars or cell array
if length (xTicks) == size (xLabels, 1),
    if (~ iscell (xLabels)) && (~ isempty (str2num (xLabels))) % labels are numeric
        bases = xTicks ./ transpose (str2num (xLabels));
        [valMin,idxMin] = min (bases);
        if valMin == max (bases), % single value
            XBase = log10 (bases (idxMin)); % XBase equivalent to Matlab figure
        end
    end
end
% Extrace bases
fig.XLim=get(ah,'XLim');
if isempty (XBase), % extrace base manually
    if log10(fig.XLim(2)-fig.XLim(1)) > 3 % round towards infinity for positiv bases
        fig.XBase=ceil(log10(abs(fig.XLim(2)-fig.XLim(1))));
    elseif log10(fig.XLim(2)-fig.XLim(1)) < -2
        fig.XBase=floor(log10(abs(fig.XLim(2)-fig.XLim(1))));
    else
        fig.XBase=1; % standard
    end
else
    fig.XBase = XBase; % use base from Figure
end

YBase   = [];
yTicks  = get (ah, 'YTick');
yLabels = get (ah, 'YTickLabel'); % array of chars or cell array
if length (yTicks) == size (yLabels, 1),
    if (~ iscell (yLabels)) && (~ isempty (str2num (yLabels))) % labels are numeric
        bases = yTicks ./ transpose (str2num (yLabels));
        [valMin,idxMin] = min (bases);
        if valMin == max (bases), % single value
            YBase = log10 (bases (idxMin)); % XBase equivalent to Matlab figure
        end
    end
end
% Extrace bases
fig.YLim=get(ah,'YLim');
if isempty (YBase), % extrace base manually
    if log10(fig.YLim(2)-fig.YLim(1)) > 3 % round towards infinity for positiv bases
        fig.YBase=ceil(log10(abs(fig.YLim(2)-fig.YLim(1))));
    elseif log10(fig.YLim(2)-fig.YLim(1)) < -2
        fig.YBase=floor(log10(abs(fig.YLim(2)-fig.YLim(1))));
    else
        fig.YBase=1; % standard
    end
else
    fig.YBase = YBase; % use base from Figure
end

% extract grids and ticks
if strcmp(get(ah,'XGrid'),'on'),        fig.XGrid='on';        else fig.XGrid='off';        end
if strcmp(get(ah,'XMinorGrid'),'on'),   fig.XMinorGrid='on';   else fig.XMinorGrid='off';   end
if strcmp(get(ah,'XMinorTick'),'on'),   fig.XMinorTick='on';   else fig.XMinorTick='off';   end
if strcmp(get(ah,'YGrid'),'on'),        fig.YGrid='on';        else fig.YGrid='off';        end
if strcmp(get(ah,'YMinorGrid'),'on'),   fig.YMinorGrid='on';   else fig.YMinorGrid='off';   end
if strcmp(get(ah,'YMinorTick'),'on'),   fig.YMinorTick='on';   else fig.YMinorTick='off';   end
fig.GridLineStyle=get(ah,'gridlinestyle'); fig.GridLineStyleTeX = convert_lineStyle(fig.GridLineStyle);
fig.MinorGridLineStyle=get(ah,'minorgridlinestyle'); fig.MinorGridLineStyleTeX = convert_lineStyle(fig.MinorGridLineStyle);
[legh,~,~,fig.LegendStr] = legend(ah); % get legend information
legend_org = legh; % copy original legend
if not(isempty(legh)), % if a legend is plotted
    fig.LegendLoc=get(legh,'Location'); % get location
    %     if strcmp(fig.LegendLoc,'none'), % no correct poition, ...
    %         fig.LegendLoc='NorthWest';   % place NorthWest
    %     end
else
    fig.LegendLoc=[];
end

fig.XLabel = get(get(ah,'XLabel'),'String'); % save xlabel
fig.YLabel = get(get(ah,'YLabel'),'String'); % save ylabel
fig.Title  = get(get(ah,'Title'),'String');
if (size(fig.XLabel,1)>1) || (size(fig.YLabel,1)>1) || (size(fig.Title,1)>1)
    h=warndlg('Multiline labels or titles are not supported at the moment. Will keep the first line only!','Warning!','modal');
    uiwait(h);
end
if size(fig.XLabel,1)>1,
    %     for ii=2:size(fig.XLabel,1), % if more than one row, convert to one row
    %         tmpstr1=[fig.XLabel(ii-1,:),'\\',fig.XLabel(ii,:)];
    %     end
    %     fig.XLabel=tmpstr1;
    fig.XLabel=fig.XLabel(1,:);
end
if size(fig.YLabel,1)>1,
    %     for ii=2:size(fig.YLabel,1), % if more than one row, convert to one row
    %         tmpstr2=[fig.YLabel(ii-1,:),'\\',fig.YLabel(ii,:)];
    %     end
    %     fig.YLabel=tmpstr2;
    fig.YLabel=fig.YLabel(1,:);
end
if size(fig.Title,1)>1,
    %     for ii=2:size(fig.Title,1), % if more than one row, convert to one row
    %         tmpstr3=[fig.Title(ii-1,:),'\\',fig.Title(ii,:)];
    %     end
    %     fig.Title=tmpstr3;
    fig.Title=fig.Title(1,:);
end

lshandles = get(ah,'Children'); % get all children of axes
lsh = double(find(handle(lshandles),'-class','graph2d.lineseries')); % find plots
lsh = [lsh; double(find(handle(lshandles),'-class','line'))]; % append 'line' plots
lsh = flipud(lsh); % flip lsh in order to get the correct order

fig.LineStyle    = get(lsh,'LineStyle'); % Linestyle
if not(iscell(fig.LineStyle)), fig.LineStyle={fig.LineStyle}; end
fig.LineStyleTeX = convert_lineStyle(fig.LineStyle);
fig.LineWidth    = get(lsh,'LineWidth'); % LineWidth
if not(iscell(fig.LineWidth)), fig.LineWidth={fig.LineWidth}; end
fig.Marker       = get(lsh,'Marker'); % Marker
if not(iscell(fig.Marker)), fig.Marker={fig.Marker}; end
fig.LineColor    = get(lsh,'Color'); % Linecolor
if not(iscell(fig.LineColor)), fig.LineColor={fig.LineColor}; end
fig.numLines=length(lsh);
fig.FaceColor = {};
for ll = 1 : fig.numLines
    fig.FaceColor {end+1} = []; % empty facecolor field, this is used
    % to check if the current line is part of a patch
end
fig.xdataLim=[inf,-inf]; % init
fig.ydataLim=[inf,-inf];
for nl=1:fig.numLines, % plot all children of axes again in linear scale
    fig.xdata{nl} = get(lsh(nl),'XData'); % get data
    fig.ydata{nl} = get(lsh(nl),'YData');
    if fig.xLog, fig.xdata{nl} = log10(fig.xdata{nl}); end % make linear
    if fig.yLog, fig.ydata{nl} = log10(fig.ydata{nl}); end % make linear
    xmin=min(fig.xdata{nl}); xmax=max(fig.xdata{nl});
    if xmin<fig.xdataLim(1), fig.xdataLim(1)=xmin; end % update if smaller
    if xmax>fig.xdataLim(2), fig.xdataLim(2)=xmax; end % update if greater
    ymin=min(fig.ydata{nl}); ymax=max(fig.ydata{nl});
    if ymin<fig.ydataLim(1), fig.ydataLim(1)=ymin; end % update if smaller
    if ymax>fig.ydataLim(2), fig.ydataLim(2)=ymax; end % update if greater
end


% Evaluate patches (Patches are treated as lines)
lsh = double (find (handle (lshandles), '-class', 'patch'));
lsh = flipud (lsh); % flip lsh in order to get the correct order
numPatches = length(lsh);

if ~ isempty (lsh), % patches have been found, append at the beginning
    
    % Linestyle
    LineStyleP = get(lsh,'LineStyle');
    if not(iscell(LineStyleP)), LineStyleP={LineStyleP}; end
    for ll = 1 : length(fig.LineStyle)
        LineStyleP{end+1, 1} = fig.LineStyle{ll};
    end
    fig.LineStyle    = LineStyleP;
    fig.LineStyleTeX = convert_lineStyle(fig.LineStyle);
    
    % Linewidth
    LineWidthP = get(lsh,'LineWidth');
    if not(iscell(LineWidthP)), LineWidthP={LineWidthP}; end
    for ll = 1 : length(fig.LineWidth)
        LineWidthP {end+1, 1} = fig.LineWidth{ll};
    end
    fig.LineWidth = LineWidthP;
    
    % Marker
    MarkerP       = get(lsh,'Marker'); % should be 'none'
    if not(iscell(MarkerP)), MarkerP={MarkerP}; end
    for ll = 1 : length(fig.Marker)
        MarkerP{end+1, 1} = fig.Marker{ll};
    end
    fig.Marker = MarkerP;
    
    
    % FaceColor as an additional field
    % test = get(lsh)
    FaceColorCData = get(lsh,'CData');
    if not(iscell(FaceColorCData)), FaceColorCData = {FaceColorCData}; end
    FaceColorP     = get(lsh,'FaceColor');
    if not(iscell(FaceColorP)), FaceColorP = {FaceColorP}; end
    [cmin, cmax] = caxis(ah);
    cmap = colormap; % get the current colormap
    m = size(cmap, 1); % number of different colors in the colormap
    for ll = 1 : length(FaceColorP)
        if isstr(FaceColorP{ll}), % test if FaceColor is given as a rgb value
            % compute color index according to Matlab help 'caxis'
            
            % test if still vector
            if length(FaceColorCData{ll}) > 1,
                % take first value
                FaceColorCData{ll} = FaceColorCData{ll}(1);
            end
            index = fix((FaceColorCData{ll}-cmin)/(cmax-cmin)*m)+1;
            if index > m, index = m; end
            FaceColorP{ll} = cmap(index,:);
        end
    end
    
    for ll = 1 : fig.numLines
        FaceColorP {end+1} = fig.FaceColor{ll};
    end
    fig.FaceColor = FaceColorP;
    if not(iscell(fig.FaceColor)), fig.FaceColor={fig.FaceColor}; end
    
    % EdgeColor -> LineColor
    EdgeColorP = get(lsh,'EdgeColor');
    if not(iscell(EdgeColorP)), EdgeColorP = {EdgeColorP}; end
    for ll = 1 : fig.numLines
        EdgeColorP{end+1} = fig.LineColor{ll};
    end
    fig.LineColor = EdgeColorP;
    if not(iscell(fig.LineColor)), fig.LineColor={fig.LineColor}; end
    
    % fig.xdataLim=[inf,-inf]; % init, not necessary, already defined
    % fig.ydataLim=[inf,-inf];
    % plot all children of axes again in linear scale
    for np = 1 : numPatches,
        xdata = get(lsh(np),'XData'); % get data
        Pxdata{np} = [xdata; xdata(1)]; % each patch must have the
        % first point as end point
        ydata = get(lsh(np),'YData');
        Pydata{np} = [ydata; ydata(1)];
        if fig.xLog, Pxdata{np} = log10(Pxdata{np}); end % make linear
        if fig.yLog, Pydata{np} = log10(Pydata{np}); end % make linear
        xmin=min(Pxdata{np}); xmax=max(Pxdata{np});
        if xmin<fig.xdataLim(1), fig.xdataLim(1)=xmin; end % update if smaller
        if xmax>fig.xdataLim(2), fig.xdataLim(2)=xmax; end % update if greater
        ymin=min(Pydata{np}); ymax=max(Pydata{np});
        if ymin<fig.ydataLim(1), fig.ydataLim(1)=ymin; end % update if smaller
        if ymax>fig.ydataLim(2), fig.ydataLim(2)=ymax; end % update if greater
    end
    
    for nl=1:fig.numLines, % plot all children of axes again in linear scale
        Pxdata{end+1} = fig.xdata{nl};
        Pydata{end+1} = fig.ydata{nl};
    end
    fig.xdata = Pxdata;
    fig.ydata = Pydata;
    fig.numLines = fig.numLines + numPatches;
    
    
end


fig.XTick=get(ah,'XTick');
fig.XTickLabel=[]; % init
if strcmp (get (ah, 'XTickLabelMode'), 'auto'),
    for tck=1:length(fig.XTick), % set all latex strings for XTickLabel
        fig.XTickLabel{tck}    = [num2str(fig.XTick(tck))];
        fig.XTickLabelTex{tck} = ['$',num2str(fig.XTick(tck)),'$'];
    end
else
    fig.XTickLabel    = convert_TickLabels2Cell (fig.XTick, get(ah,'XTickLabel'));
    if ~ text_contains_dollarsymbol (fig.XTickLabel)
        fig.XTickLabelTex = convert_text_to_latex (fig.XTickLabel);
    else
        fig.XTickLabelTex = fig.XTickLabel; % already TeX format
    end
end

if fig.xLog,
    XTickNew = log10(fig.XTick);
    XLimNew  = log10(fig.XLim);
    if XLimNew(2) - XLimNew(1) < 1, % too small for logarithmic scale
        XTickNew = [floor(fig.xdataLim(1)):1:ceil(fig.xdataLim(2))];
        XLimNew  = [floor(fig.xdataLim(1)),ceil(fig.xdataLim(2))];
    end
    if (length(XTickNew) == length(fig.XTick)) && ... % maintain old labels
            (strcmp (get (ah, 'XTickLabelMode'), 'manual'))
        % do nothing
    else
        for tck=1:length(fig.XTick), % update latex strings for XTickLabel
            fig.XTickLabel{tck}    = ['10^',num2str(XTickNew(tck))];
            fig.XTickLabelTex{tck} = ['$10^{',num2str(XTickNew(tck)),'}$'];
        end
    end
    fig.XTick = XTickNew;
    fig.XLim  = XLimNew;
end



fig.YTick=get(ah,'YTick');
fig.YTickLabel=[]; % init
if strcmp (get (ah, 'YTickLabelMode'), 'auto'),
    for tck=1:length(fig.YTick), % set all latex strings for XTickLabel
        fig.YTickLabel{tck}    = [num2str(fig.YTick(tck))];
        fig.YTickLabelTex{tck} = ['$',num2str(fig.YTick(tck)),'$'];
    end
else
    fig.YTickLabel    = convert_TickLabels2Cell (fig.YTick, get(ah,'YTickLabel'));
    if ~ text_contains_dollarsymbol (fig.YTickLabel)
        fig.YTickLabelTex = convert_text_to_latex (fig.YTickLabel);
    else
        fig.YTickLabelTex = fig.YTickLabel; % already TeX format
    end
end

if fig.yLog,
    YTickNew = log10(fig.YTick);
    YLimNew  = log10(fig.YLim);
    if YLimNew(2) - YLimNew(1) < 1, % too small for logarithmic scale
        YTickNew = [floor(fig.ydataLim(1)):1:ceil(fig.ydataLim(2))];
        YLimNew  = [floor(fig.ydataLim(1)),ceil(fig.ydataLim(2))];
    end
    if (length(YTickNew) == length(fig.YTick)) && ... % maintain old labels
            (strcmp (get (ah, 'YTickLabelMode'), 'manual'))
        % do nothing
    else
        for tck=1:length(fig.YTick), % update latex strings for XTickLabel
            fig.YTickLabel{tck}    = ['10^']; % no exponent, dummy, just for visualization
            fig.YTickLabelTex{tck} = ['$10^{',num2str(YTickNew(tck)),'}$'];
        end
    end
    fig.YTick = YTickNew;
    fig.YLim  = YLimNew;
end



% Test for too extrem X-values in linear case
fig.XLabelPost = '';
if not(fig.xLog) && (fig.XBase>2 || fig.XBase<-2), % if not logarithmic
    % change the following variables:
    % xdata, xdataLim, XTick, XTickLabel, XTickLabelTex, XLim
    for nl=1:fig.numLines, % plot all children of axes again in linear scale
        fig.xdata{nl}=fig.xdata{nl}/(10^fig.XBase);
    end
    fig.xdataLim=fig.xdataLim/(10^fig.XBase);
    fig.XTick=fig.XTick/(10^fig.XBase);
    fig.XLim=fig.XLim/(10^fig.XBase);
    if strcmp (get (ah, 'XTickLabelMode'), 'auto'),
        for tck=1:length(fig.XTick), % set all latex strings for XTickLabel
            fig.XTickLabel{tck}    = [num2str(fig.XTick(tck))];%,'e',num2str(fig.XBase)];
            fig.XTickLabelTex{tck} = ['$',num2str(fig.XTick(tck)),'$'];
        end
        fig.XLabelPost = [' $ \cdot 10^{',num2str(fig.XBase),'}$'];
    end
end

% Test for too extrem Y-values in linear case
fig.YLabelPost = '';
if not(fig.yLog) && (fig.YBase>2 || fig.YBase<-2), % if not logarithmic
    % change the following variables:
    % ydata, ydataLim, YTick, YTickLabel, YTickLabelTex, YLim
    for nl=1:fig.numLines, % plot all children of axes again in linear scale
        fig.ydata{nl}=fig.ydata{nl}/(10^fig.YBase);
    end
    fig.ydataLim=fig.ydataLim/(10^fig.YBase);
    fig.YTick=fig.YTick/(10^fig.YBase);
    fig.YLim=fig.YLim/(10^fig.YBase);
    if strcmp (get (ah, 'YTickLabelMode'), 'auto'),
        for tck=1:length(fig.YTick), % set all latex strings for YTickLabel
            fig.YTickLabel{tck}    = [num2str(fig.YTick(tck))];%,'e',num2str(fig.YBase)];
            fig.YTickLabelTex{tck} = ['$',num2str(fig.YTick(tck)),'$'];
        end
        fig.YLabelPost = [' $ \cdot 10^{',num2str(fig.YBase),'}$'];
    end
end

% extract annotations and convert the x- and y- values
% Textarrows
fig.TextarrowsStr=[]; %init
allhandles = findall(fh);
anhar = double(find(handle(allhandles),'-class','scribe.textarrow'));
anhar = anhar(1:end/2); % take only the first half
for an= 1:length(anhar),
    text_ta = get(anhar(an),'String');
    if ischar(text_ta) % for a specific reason that I can't figure out,
        % it is a char array here, -> convert
        for jj = 1: size(text_ta,1), % all rows
            text_ta_out{jj,1} = text_ta(jj,:);
        end
        fig.TextarrowsStr{an}   = text_ta_out; % copy
    else
        fig.TextarrowsStr{an}   = text_ta;
    end
    fig.TextarrowsX(an,:)   = get(anhar(an),'X'); % normalized figure units
    fig.TextarrowsY(an,:)   = get(anhar(an),'Y');
    fig.TextarrowsRot(an)   = 180 + 180/pi ...
        * atan2(fig.TextarrowsY(an,2)-fig.TextarrowsY(an,1), ...
        fig.TextarrowsX(an,2)-fig.TextarrowsX(an,1));
    fig.TextarrowsColor{an} = get(anhar(an),'color');
    fig.TextarrowsLineStyle{an} = get(anhar(an),'linestyle');
    fig.TextarrowsHorizontalAlignment{an} = get(anhar(an),'HorizontalAlignment');
end
if ~isempty(anhar)
    % convert relative figure positions to x-and y-values
    set(ah,'Units','normalized')
    pos=get(ah,'Position');
    XY1=pos(1:2); Width=pos(3:4);
    fig.TextarrowsX=(fig.TextarrowsX-XY1(1))./Width(1); % normalized axes units
    fig.TextarrowsY=(fig.TextarrowsY-XY1(2))./Width(2);
    fig.TextarrowsX=fig.TextarrowsX.*(fig.XLim(2)-fig.XLim(1))+fig.XLim(1); % matched x-units
    fig.TextarrowsY=fig.TextarrowsY.*(fig.YLim(2)-fig.YLim(1))+fig.YLim(1); % matched y-units
end
% TextBoxes
fig.TextBoxesStr=[]; %init
anhtb = double(find(handle(allhandles),'-class','scribe.textbox'));
anhtb = anhtb(1:end/2); % take only the first half
for an= 1:length(anhtb),
    text_tb = get(anhtb(an),'String');
    if ischar(text_tb) % for a specific reason that I can't figure out,
        % it is a char array here, -> convert
        for jj = 1: size(text_tb,1), % all rows
            text_tb_out{jj,1} = text_tb(jj,:);
        end
        fig.TextBoxesStr{an}   = text_tb_out; % copy
    else
        fig.TextBoxesStr{an}   = text_tb;
    end
    fig.TextBoxesPosition(an,:)   = get(anhtb(an),'Position'); % normalized figure units
    fig.TextBoxesColor{an} = get(anhtb(an),'Color');
    fig.TextBoxesLineStyle{an} = get(anhtb(an),'LineStyle');
    fig.TextBoxesHorizontalAlignment{an} = get(anhtb(an),'HorizontalAlignment');
end
if ~isempty(anhtb)
    % convert relative figure positions to x-and y-values
    set(ah,'Units','normalized')
    pos=get(ah,'Position');
    XY1=pos(1:2); Width=pos(3:4);
    fig.TextBoxesX = (fig.TextBoxesPosition(:,1)-XY1(1))./Width(1); % normalized axes units
    fig.TextBoxesY = (fig.TextBoxesPosition(:,2)-XY1(2))./Width(2);
    fig.TextBoxesX = fig.TextBoxesX.*(fig.XLim(2)-fig.XLim(1))+fig.XLim(1); % matched x-units
    fig.TextBoxesY = fig.TextBoxesY.*(fig.YLim(2)-fig.YLim(1))+fig.YLim(1); % matched y-units
end

% extract box
if strcmp(get(ah,'Box'),'on'), fig.Box='on'; else fig.Box='off'; end

% if none of the strings contains a '$', set converting on by default
if  ~ text_contains_dollarsymbol (fig.LegendStr)
    fig.LegendStr = convert_text_to_latex (fig.LegendStr);
end
if  ~ text_contains_dollarsymbol (fig.XLabel)
    fig.XLabel = convert_text_to_latex (fig.XLabel);
end
if  ~ text_contains_dollarsymbol (fig.YLabel)
    fig.YLabel = convert_text_to_latex (fig.YLabel);
end
if  ~ text_contains_dollarsymbol (fig.Title)
    fig.Title = convert_text_to_latex (fig.Title);
end
for kk = 1 : length (fig.TextarrowsStr)
    if ~text_contains_dollarsymbol(fig.TextarrowsStr{kk})
        fig.TextarrowsStr{kk} = convert_text_to_latex (fig.TextarrowsStr{kk});
    end
end
for kk = 1 : length (fig.TextBoxesStr)
    if ~text_contains_dollarsymbol(fig.TextBoxesStr{kk})
        fig.TextBoxesStr{kk} = convert_text_to_latex (fig.TextBoxesStr{kk});
    end
end

%% Extract Data of second yaxis
fig.YMinorTick2='off'; % init
if ~isempty(ahSecond), % process second y-axis
    % all variables of the second axis are same as first with an extra "2":
    
    % Test if Scales are log
    fig.xLog2=fig.xLog; % x-axis is equal
    if strcmp(get(ahSecond,'YScale'),'log'), fig.yLog2=1; else fig.yLog2=0; end
    
    % Extrace bases
    fig.XLim2=fig.XLim;
    fig.XBase2=fig.XBase;
    YBase2  = [];
    yTicks  = get (ahSecond, 'YTick');
    yLabels = get (ahSecond, 'YTickLabel'); % array of chars or cell array
    if length (yTicks) == size (yLabels, 1),
        if (~ iscell (yLabels)) && (~ isempty (str2num (yLabels))) % labels are numeric
            bases = yTicks ./ transpose (str2num (yLabels));
            [valMin,idxMin] = min (bases);
            if valMin == max (bases), % single value
                YBase2 = log10 (bases (idxMin)); % YBase equivalent to Matlab figure
            end
        end
    end
    fig.YLim2=get(ahSecond,'YLim');
    if isempty (YBase2), % extrace base manually
        if log10(fig.YLim2(2)-fig.YLim2(1)) > 3 % round towards infinity for positiv bases
            fig.YBase2=ceil(log10(abs(fig.YLim2(2)-fig.YLim2(1))));
        elseif log10(fig.YLim2(2)-fig.YLim2(1)) < -2
            fig.YBase2=floor(log10(abs(fig.YLim2(2)-fig.YLim2(1))));
        else
            fig.YBase2=1; % standard
        end
    else
        fig.YBase2 = YBase2; % use base from Figure
    end
    
    % extract grids and ticks
    fig.XGrid2='off';
    fig.XMinorGrid2='off';
    fig.XMinorTick2='off';
    fig.XTick2='off';
    fig.YGrid2='off';
    fig.YMinorGrid2='off';
    if strcmp(get(ahSecond,'YMinorTick'),'on'),   fig.YMinorTick2='on';   else fig.YMinorTick2='off';   end
    fig.LegendLoc2=[]; % no second legend, all information is in first one
    fig.XLabel2 = '';
    fig.YLabel2 = get(get(ahSecond,'YLabel'),'String'); % save ylabel
    if  ~ text_contains_dollarsymbol (fig.YLabel2)
        fig.YLabel2 = convert_text_to_latex (fig.YLabel2);
    end
    
    % extrace data of second axis
    lshandles2 = get(ahSecond,'Children'); % get all children of axes
    lsh2 = double(find(handle(lshandles2),'-class','graph2d.lineseries'));
    lsh2 = [lsh2; double(find(handle(lshandles2),'-class','line'))]; % append 'line' plots
    lsh2 = flipud(lsh2); % flip lsh in order to get the correct order
    fig.LineStyle2    = get(lsh2,'LineStyle'); % Linestyle
    if not(iscell(fig.LineStyle2)), fig.LineStyle2={fig.LineStyle2}; end
    fig.LineStyleTeX2 = convert_lineStyle(fig.LineStyle2);
    fig.LineWidth2    = get(lsh2,'LineWidth'); % LineWidth
    if not(iscell(fig.LineWidth2)), fig.LineWidth2={fig.LineWidth2}; end
    fig.Marker2       = get(lsh2,'Marker'); % Marker
    if not(iscell(fig.Marker2)), fig.Marker2={fig.Marker2}; end
    fig.LineColor2    = get(lsh2,'Color'); % Linecolor
    if not(iscell(fig.LineColor2)), fig.LineColor2={fig.LineColor2}; end
    fig.numLines2=length(lsh2);
    for ll = 1 : fig.numLines2
        fig.FaceColor {end+1} = []; % empty facecolor field,
    end
    fig.xdataLim2=fig.xdataLim;
    fig.ydataLim2=[inf,-inf]; % init
    for nl=1:fig.numLines2, % plot all children of axes again in linear scale
        fig.xdata2{nl} = get(lsh2(nl),'XData'); % get data
        fig.ydata2{nl} = get(lsh2(nl),'YData');
        if fig.xLog2, fig.xdata2{nl} = log10(fig.xdata2{nl}); end % make linear
        if fig.yLog2, fig.ydata2{nl} = log10(fig.ydata2{nl}); end % make linear
        ymin=min(fig.ydata2{nl}); ymax=max(fig.ydata2{nl});
        if ymin<fig.ydataLim2(1), fig.ydataLim2(1)=ymin; end % update if smaller
        if ymax>fig.ydataLim2(2), fig.ydataLim2(2)=ymax; end % update if greater
    end
    
    fig.YTick2=get(ahSecond,'YTick');
    fig.YTickLabel2=[]; % init
    if strcmp (get (ahSecond, 'YTickLabelMode'), 'auto'),
        for tck=1:length(fig.YTick2), % set all latex strings for XTickLabel
            fig.YTickLabel2{tck}    = [num2str(fig.YTick2(tck))];
            fig.YTickLabelTex2{tck} = ['$',num2str(fig.YTick2(tck)),'$'];
        end
    else
        fig.YTickLabel2    = convert_TickLabels2Cell (fig.YTick2, get(ahSecond,'YTickLabel'));
        if ~ text_contains_dollarsymbol (fig.YTickLabel2)
            fig.YTickLabelTex2 = convert_text_to_latex (fig.YTickLabel2);
        else
            fig.YTickLabelTex2 = fig.YTickLabel2; % already TeX format
        end
    end
    if fig.yLog2,
        fig.YTick2=log10(fig.YTick2);
        fig.YLim2=log10(fig.YLim2);
        if fig.YLim2(2)-fig.YLim2(1)<1, % too small for logarithmic scale
            fig.YTick2=[floor(fig.ydataLim2(1)):1:ceil(fig.ydataLim2(2))];
            fig.YLim2=[floor(fig.ydataLim2(1)),ceil(fig.ydataLim2(2))];
        end
        for tck=1:length(fig.YTick2), % set all latex strings for XTickLabel
            fig.YTickLabel2{tck}    = ['10^']; % no exponent, dummy, just for visualization
            fig.YTickLabelTex2{tck} = ['$10^{',num2str(fig.YTick2(tck)),'}$'];
        end
    end
    
    % Test for too extrem X-values in linear case
    if not(fig.xLog) && (fig.XBase>2 || fig.XBase<-1), % if not logarithmic
        % change xdata2:
        for nl=1:fig.numLines2, % plot all children of axes again in linear scale
            fig.xdata2{nl}=fig.xdata2{nl}/(10^fig.XBase);
        end
    end
    
    % Test for too extrem Y-values in linear case
    fig.YLabelPost2 = '';
    if not(fig.yLog2) && (fig.YBase2>2 || fig.YBase2<-2), % if not logarithmic
        % change the following variables:
        % ydata, ydataLim, YTick, YTickLabel, YTickLabelTex, YLim
        for nl=1:fig.numLines2, % plot all children of axes again in linear scale
            fig.ydata2{nl}=fig.ydata2{nl}/(10^fig.YBase2);
        end
        fig.ydataLim2=fig.ydataLim2/(10^fig.YBase2);
        fig.YTick2=fig.YTick2/(10^fig.YBase2);
        fig.YLim2=fig.YLim2/(10^fig.YBase2);
        if strcmp (get (ahSecond, 'YTickLabelMode'), 'auto'),
            for tck=1:length(fig.YTick2), % set all latex strings for YTickLabel
                fig.YTickLabel2{tck}    = [num2str(fig.YTick2(tck))];%,'e',num2str(fig.YBase)];
                fig.YTickLabelTex2{tck} = ['$',num2str(fig.YTick2(tck)),'$'];
            end
            fig.YLabelPost2 = [' $ \cdot 10^{',num2str(fig.YBase2),'}$'];
        end
    end
    
    % convert y-data, and y-ticks
    for nl=1:fig.numLines2, % convert all y-data in order to match with the first axes
        fig.ydata2{nl}=((fig.ydata2{nl}-fig.YLim2(1))/(fig.YLim2(2)-fig.YLim2(1))*...
            (fig.YLim(2)-fig.YLim(1)))+fig.YLim(1);
    end
    fig.YTick2Tex=((fig.YTick2-fig.YLim2(1))/(fig.YLim2(2)-fig.YLim2(1))*...
        (fig.YLim(2)-fig.YLim(1)))+fig.YLim(1);
    yAxisCenversionFac=(fig.YLim(2)-fig.YLim(1))/(fig.YLim2(2)-fig.YLim2(1));
    
    % append data
    for nl=1:fig.numLines2, % append data
        fig.xdata{fig.numLines+nl}=fig.xdata2{nl};
        fig.ydata{fig.numLines+nl}=fig.ydata2{nl};
        fig.LineColor{fig.numLines+nl}=fig.LineColor2{nl};
        fig.LineStyle{fig.numLines+nl}=fig.LineStyle2{nl};
        fig.LineWidth{fig.numLines+nl}=fig.LineWidth2{nl};
        fig.Marker{fig.numLines+nl}=fig.Marker2{nl};
    end
    fig.numLines=fig.numLines+fig.numLines2; % update number of lines
end


%% PLOT EXTRACTED DATA IN NEW FIGURE

% check for docked DefaultFigureWindowStyle
if strcmp (get(0,'DefaultFigureWindowStyle'), 'docked')
    warning(['fig2TexPS won''t work with docked figures correctly. ', ...
        'Please switch of this behavior by typing: ', ...
        'set(0,''DefaultFigureWindowStyle'',''normal''). ', ...
        'The docked behavior can be set afterwards back to: ',...
        'set(0,''DefaultFigureWindowStyle'',''docked'')'])
end

% create new figure and new axes that will have linear axes:
fig.figh=figure('Name','Figure outline that will be converted to LaTeX',...
    'HitTest','off','MenuBar','none',...
    'NumberTitle','off','Resize','off');
fig.axh=axes;

for nl=1:fig.numLines,
    plot(fig.axh,fig.xdata{nl},fig.ydata{nl},...
        'Color',fig.LineColor{nl},...
        'LineStyle',fig.LineStyle{nl},...
        'LineWidth',fig.LineWidth{nl},...
        'Marker',fig.Marker{nl});
    hold on;
end
set(fig.axh, ...
    'XGrid',fig.XGrid,...
    'YGrid',fig.YGrid,...
    'Box',fig.Box,...
    'gridlinestyle',fig.GridLineStyle,...
    'minorgridlinestyle',fig.MinorGridLineStyle,...
    'Visible',fig.AxesVisible);
if not(fig.xLog), % if not logarithmic, set minor grid
    set(fig.axh,'XMinorGrid',fig.XMinorGrid);
    set(fig.axh,'XMinorTick',fig.XMinorTick);
end
if not(fig.yLog),
    set(fig.axh,'YMinorGrid',fig.YMinorGrid);
    set(fig.axh,'YMinorTick',fig.YMinorTick);
end
set(fig.axh,'XTick',fig.XTick,'XLim',fig.XLim);
set(fig.axh,'XTickLabel',fig.XTickLabel)
set(fig.axh,'YTick',fig.YTick,'YLim',fig.YLim);
set(fig.axh,'YTickLabel',fig.YTickLabel)
set(fig.axh,'FontSize',globals.labelsep)
xlabel(fig.XLabel,'FontSize',globals.labelsep,'interpreter','none');
ylabel(fig.YLabel,'FontSize',globals.labelsep,'interpreter','none');
title(fig.Title,'FontSize',globals.labelsep,'interpreter','none');
if not(isempty(fig.LegendLoc)), % a legend is replotted
    legend(fig.axh,fig.LegendStr,'Location',fig.LegendLoc,'Interpreter','none')
end
if ~isempty(ahSecond), % process second y-axis
    fig.axh2=axes('Units',get(fig.axh,'Units'),...
        'Position',get(fig.axh,'Position'),...
        'Parent',get(fig.axh,'Parent'),...
        'XAxisLocation','bottom','YAxisLocation','right',...
        'Color','none','XColor','k','XTick',[],'YTickLabel',fig.YTickLabel2,...
        'YColor','k','YLim',fig.YLim2,'YTick',fig.YTick2,'YMinorTick',fig.YMinorTick2);
    ylabel(fig.YLabel2,'FontSize',globals.labelsep,'interpreter','none');
end

movegui(fig.figh,'center'); % shift figure to screen center



%% FIGURE GEOMETRY SETTINGS
% the legend is turned off inside this function if it is outside the graph
% otherwise it is showed in the figure
legend(fig.axh,'off'); % turn off the legend
[fig,globals] = update_figure_geometry(fig,globals); % subfunction GUI

if exit_flag,
    warning(s); % turn on specified warnings again
    return,
end


fp=get(fig.figh,'Position');

set(fig.axh,'Units','normalized'); ap=get(fig.axh,'Position'); %[left bottom width height]
ti=get(fig.axh,'TightInset'); % [left bottom right top]

if ~isempty(ahSecond),
    set(fig.axh2,'Units','normalized');
    ti2=get(fig.axh2,'TightInset'); % [left bottom right top]
    ti(3)=ti2(3); % update extra space on the right side caused by second y-label
end

if not(isempty(fig.LegendLoc)) && ~isempty(regexp(fig.LegendLoc,'Outside', 'once')), % a legend is replotted if the legend is outside
    legend(fig.axh,fig.LegendStr,'Location',fig.LegendLoc,'interpreter','none')
end


% test if axes should be centered, if yes, adapt either left ti(1) or right
% ti(1) to the greater one
if globals.centerAxes==1,
    if ti(3)<ti(1), ti(3)=ti(1); else ti(1)=ti(3); end
end
% apPlusTight=[ap(1)-ti(1),ap(2)-ti(2),ap(3)+ti(1)+ti(3),ap(4)+ti(2)+ti(4)];

DeltaXLim=fig.XLim(2)-fig.XLim(1);
xunit = 1/DeltaXLim;
DeltaYLim=fig.YLim(2)-fig.YLim(1);
yunit = globals.heightFactor/DeltaYLim;
% Here xunit and yunit are set in order that the width matches to the axes
% If the width of the figure shoud match to the figure, both xunit and
% yunit must be scaled down
if strcmp(globals.widthType,'figure'),
    % ap(3)       % normalized width of axes
    % ti(1)+ti(3) % normalized width of extra space
    totWidthNorm=ap(3)+ti(1)+ti(3);
    xunit=xunit*ap(3)/totWidthNorm; % scale
    yunit=yunit*ap(3)/totWidthNorm; % scale
%     StrFile{end+1}=['\providelength{\plotwidth}           \setlength{\plotwidth}{',num2str(globals.plotwidth),'cm}% width of the total figure']; % append next line
else
%     StrFile{end+1}=['\providelength{\plotwidth}           \setlength{\plotwidth}{',num2str(globals.plotwidth),'cm}% width of the axes only']; % append next line
end
% create the length, with options to re-write it.
if isnumeric(globals.plotwidth)
    newWidth = [num2str(globals.plotwidth),'cm'];
else
    newWidth = globals.plotwidth;
end
StrFile{end+1}='\ifx \plotwidth \undefined% if this is the first figure';
StrFile{end+1}=['    \providelength{\plotwidth} \setlength{\plotwidth}{',newWidth,'}%'];
StrFile{end+1}='\else% this figure is not the only one in the file';
StrFile{end+1}='  \ifdim \plotwidth=0pt%';
StrFile{end+1}=['    \setlength{\plotwidth}{',newWidth,'}%'];
StrFile{end+1}='  \fi%';
StrFile{end+1}='\fi%';


if ischar(globals.LineWidth)
    StrFile{end+1}=['\providelength{\LineWidth}           \setlength{\LineWidth}{',globals.LineWidth,'}%']; % append next line
else
    StrFile{end+1}=['\providelength{\LineWidth}           \setlength{\LineWidth}{',num2str(globals.LineWidth),'pt}%']; % append next line
end
if numPatches>0,
    if ischar(globals.EdgeLineWidth)
        StrFile{end+1}=['\providelength{\EdgeLineWidth}       \setlength{\EdgeLineWidth}{',globals.EdgeLineWidth,'}%']; % append next line
    else
        StrFile{end+1}=['\providelength{\EdgeLineWidth}       \setlength{\EdgeLineWidth}{',num2str(globals.EdgeLineWidth),'pt}%']; % append next line
    end
end


if ischar(globals.MarkerSize)
    StrFile{end+1}=['\providelength{\MarkerSize}          \setlength{\MarkerSize}{',globals.MarkerSize,'}%']; % append next line
else
    StrFile{end+1}=['\providelength{\MarkerSize}          \setlength{\MarkerSize}{',num2str(globals.MarkerSize),'pt}%']; % append next line
end
StrFile{end+1}=['\newrgbcolor{GridColor}{',num2str(globals.GridGreyFac),' ',num2str(globals.GridGreyFac),' ',num2str(globals.GridGreyFac),'}%'];
StrFile{end+1}='%';
StrFile{end+1}='% Begin Figure:-------------------------------------------';
tmpstr = '\psset{xunit=#1#\plotwidth,yunit=#2#\plotwidth}%';
StrFile{end+1}=strfill(tmpstr,{xunit,yunit});

if ~isempty(globals.tightInset_cm) && strcmp(globals.widthType,'axes'),
    % Here, the fixed borders are set
    figx1 = fig.XLim(1)-globals.tightInset_cm(1)/globals.plotwidth*DeltaXLim;
    figy1 = fig.YLim(1)-globals.tightInset_cm(2)/globals.plotwidth*DeltaYLim/globals.heightFactor;
    figx2 = fig.XLim(1)+DeltaXLim+globals.tightInset_cm(3)/globals.plotwidth*DeltaXLim;
    figy2 = fig.YLim(1)+DeltaYLim+globals.tightInset_cm(4)/globals.plotwidth*DeltaYLim/globals.heightFactor;
else
    if ~isempty(globals.tightInset_cm),
        % Here, conflict because of fixed figure size
        warning('The extra space setting will work correctly only when the width of the axes is specified!')
    end
    figx1 = fig.XLim(1)-ti(1)/ap(3)*DeltaXLim;
    figy1 = fig.YLim(1)-ti(2)/ap(4)*DeltaYLim;
    figx2 = fig.XLim(1)+(ap(3)+ti(3))/ap(3)*DeltaXLim;
    figy2 = fig.YLim(1)+(ap(4)+ti(4))/ap(4)*DeltaYLim;
end


tmpstr = '\begin{pspicture}(#1#,#2#)(#3#,#4#)%';
StrFile{end+1}=strfill(tmpstr,{figx1,figy1,figx2,figy2});

%% DRAW GRID AND AXES
% draw box for visualization
StrFile{end+1}='';
StrFile{end+1}='% Draw bounding box for test aspects: ----';
if globals.BoundingBox
    tmpstr = '\psframe(#1#,#2#)(#3#,#4#)';
else
    tmpstr = '% \psframe(#1#,#2#)(#3#,#4#)'; % with comments
end
StrFile{end+1}=strfill(tmpstr,{figx1,figy1,figx2,figy2});
tmpstr = '% Total width:  #1# cm';
TotalWidth_cm = globals.plotwidth*xunit*(figx2-figx1);
StrFile{end+1}=strfill(tmpstr,{TotalWidth_cm});
tmpstr = '% Total height: #1# cm';
TotalHeight_cm = globals.plotwidth*yunit*(figy2-figy1);
StrFile{end+1}=strfill(tmpstr,{TotalHeight_cm});

TickMajorLengthX=globals.TickMajorLengthRel*DeltaXLim; % length of the grid lines
TickMinorLengthX=globals.TickMinorLengthRel*DeltaXLim;
TickMajorLengthY=globals.TickMajorLengthRel*DeltaYLim/globals.heightFactor;
TickMinorLengthY=globals.TickMinorLengthRel*DeltaYLim/globals.heightFactor;

if strcmp(fig.AxesVisible,'on'), % draw only in this case
    
    % draw grid
    GLSTex=fig.GridLineStyleTeX{1}; % convert to char
    if strcmp(fig.XGrid,'on') || strcmp(fig.YGrid,'on') % At least one Grid
        StrFile{end+1}='';
        StrFile{end+1}='% Draw Grid: ----';
        if ischar(globals.MinorGridLineDotSep)
            StrFile={['\providelength{\MinorGridLineDotSep} \setlength{\MinorGridLineDotSep}{',globals.MinorGridLineDotSep,'}%'],StrFile{1:end}};
        else
            StrFile={['\providelength{\MinorGridLineDotSep} \setlength{\MinorGridLineDotSep}{',num2str(globals.MinorGridLineDotSep),'pt}%'],StrFile{1:end}};
        end
        if ischar(globals.MinorGridLineWidth)
            StrFile={['\providelength{\MinorGridLineWidth}  \setlength{\MinorGridLineWidth}{',globals.MinorGridLineWidth,'}%'],StrFile{1:end}};
        else
            StrFile={['\providelength{\MinorGridLineWidth}  \setlength{\MinorGridLineWidth}{',num2str(globals.MinorGridLineWidth),'pt}%'],StrFile{1:end}};
        end
        if ischar(globals.GridLineDotSep)
            StrFile={['\providelength{\GridLineDotSep}      \setlength{\GridLineDotSep}{',globals.GridLineDotSep,'}%'],StrFile{1:end}};
        else
            StrFile={['\providelength{\GridLineDotSep}      \setlength{\GridLineDotSep}{',num2str(globals.GridLineDotSep),'pt}%'],StrFile{1:end}};
        end
        if ischar(globals.GridLineWidth)
            StrFile={['\providelength{\GridLineWidth}       \setlength{\GridLineWidth}{',globals.GridLineWidth,'}%'],StrFile{1:end}};
        else
            StrFile={['\providelength{\GridLineWidth}       \setlength{\GridLineWidth}{',num2str(globals.GridLineWidth),'pt}%'],StrFile{1:end}};
        end
        if strcmp(fig.XGrid,'on') % xgrid
            StrFile{end+1}='% x-Grid:';
            for tck=1:length(fig.XTick), % x=const
                StrFile{end+1}=['\psline[linestyle=',GLSTex,',dotsep=\GridLineDotSep,linewidth=\GridLineWidth,linecolor=GridColor](',num2str(fig.XTick(tck),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')(',num2str(fig.XTick(tck),globals.FloatingPointFormat),',',num2str(fig.YLim(2),globals.FloatingPointFormat),')'];
            end
        end
        if strcmp(fig.YGrid,'on') % ygrid
            StrFile{end+1}='% y-Grid:';
            for tck=1:length(fig.YTick), % y=const
                StrFile{end+1}=['\psline[linestyle=',GLSTex,',dotsep=\GridLineDotSep,linewidth=\GridLineWidth,linecolor=GridColor](',num2str(fig.XLim(1),globals.FloatingPointFormat),',',num2str(fig.YTick(tck),globals.FloatingPointFormat),')(',num2str(fig.XLim(2),globals.FloatingPointFormat),',',num2str(fig.YTick(tck),globals.FloatingPointFormat),')'];
            end
        end
    end
    
    % draw minor grid
    subPoints_log=log10([2 3 4 5 6 7 8 9]);
    MGLSTex=fig.MinorGridLineStyleTeX{1}; % convert to char
    if strcmp(fig.XMinorGrid,'on') || strcmp(fig.YMinorGrid,'on') % At least one Grid
        StrFile{end+1}='';
        StrFile{end+1}='% Draw Minor Grid: ----';
        if strcmp(fig.XMinorGrid,'on') && fig.xLog==1 % XMinorGrid and log
            if length(fig.XTick)>1, % at least two points
                if fig.XTick(2)-fig.XTick(1) == 1 % one decade per point
                    StrFile{end+1}='% x-Minor Grid:';
                    for tck=1:length(fig.XTick)-1, % x=const
                        for jj=1:length(subPoints_log),
                            if fig.XTick(tck)+subPoints_log(jj) <= fig.XLim(2),
                                StrFile{end+1}=['\psline[linestyle=',MGLSTex,',dotsep=\MinorGridLineDotSep,linewidth=\MinorGridLineWidth,linecolor=GridColor](',num2str(fig.XTick(tck)+subPoints_log(jj),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')(',num2str(fig.XTick(tck)+subPoints_log(jj),globals.FloatingPointFormat),',',num2str(fig.YLim(2),globals.FloatingPointFormat),')'];
                            end
                        end
                    end
                end
            end
        end
        if strcmp(fig.YMinorGrid,'on') && fig.yLog==1 % YMinorGrid and log
            if length(fig.YTick)>1, % at least two points
                if fig.YTick(2)-fig.YTick(1) == 1 % one decade per point
                    StrFile{end+1}='% y-Minor Grid:';
                    for tck=1:length(fig.YTick)-1, % x=const
                        for jj=1:length(subPoints_log),
                            if fig.YTick(tck)+subPoints_log(jj) <= fig.YLim(2),
                                StrFile{end+1}=['\psline[linestyle=',MGLSTex,',dotsep=\MinorGridLineDotSep,linewidth=\MinorGridLineWidth,linecolor=GridColor](',num2str(fig.XLim(1),globals.FloatingPointFormat),',',num2str(fig.YTick(tck)+subPoints_log(jj),globals.FloatingPointFormat),')(',num2str(fig.XLim(2),globals.FloatingPointFormat),',',num2str(fig.YTick(tck)+subPoints_log(jj),globals.FloatingPointFormat),')'];
                            end
                        end
                    end
                end
            end
        end
    end
    
    % draw ticks
    StrFile{end+1}='';
    StrFile{end+1}='% Draw Ticks: ----';
    StrFile{end+1}='% x-Ticks:';
    for tck=1:length(fig.XTick), % x=const
        StrFile{end+1}=['\psline[linewidth=\AxesLineWidth,linecolor=GridColor](',num2str(fig.XTick(tck),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')(',num2str(fig.XTick(tck),globals.FloatingPointFormat),',',num2str(fig.YLim(1)+TickMajorLengthY,globals.FloatingPointFormat),')'];
    end
    StrFile{end+1}='% y-Ticks:';
    for tck=1:length(fig.YTick), % y=const
        StrFile{end+1}=['\psline[linewidth=\AxesLineWidth,linecolor=GridColor](',num2str(fig.XLim(1),globals.FloatingPointFormat),',',num2str(fig.YTick(tck),globals.FloatingPointFormat),')(',num2str(fig.XLim(1)+TickMajorLengthX,globals.FloatingPointFormat),',',num2str(fig.YTick(tck),globals.FloatingPointFormat),')'];
    end
    if ~isempty(ahSecond),
        StrFile{end+1}='% y-Ticks Second Axis:';
        for tck=1:length(fig.YTick2), % y=const
            StrFile{end+1}=['\psline[linewidth=\AxesLineWidth,linecolor=GridColor](',num2str(fig.XLim(2)-TickMajorLengthX,globals.FloatingPointFormat),',',num2str(fig.YTick2Tex(tck),globals.FloatingPointFormat),')(',num2str(fig.XLim(2),globals.FloatingPointFormat),',',num2str(fig.YTick2Tex(tck),globals.FloatingPointFormat),')'];
        end
    end
    
    % draw minor ticks
    if strcmp(fig.XMinorTick,'on') || strcmp(fig.YMinorTick,'on') || strcmp(fig.YMinorTick2,'on') % At least one Tick
        StrFile{end+1}='';
        StrFile{end+1}='% Draw Minor Ticks: ----';
        if strcmp(fig.XMinorTick,'on') && fig.xLog==1 % XMinorTick and log
            if length(fig.XTick)>1, % at least two points
                if fig.XTick(2)-fig.XTick(1) == 1 % one decade per point
                    StrFile{end+1}='% x-Minor Ticks:';
                    for tck=1:length(fig.XTick)-1, % x=const
                        for jj=1:length(subPoints_log),
                            if fig.XTick(tck)+subPoints_log(jj) <= fig.XLim(2),
                                StrFile{end+1}=['\psline[linewidth=\AxesLineWidth,linecolor=GridColor](',num2str(fig.XTick(tck)+subPoints_log(jj),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')(',num2str(fig.XTick(tck)+subPoints_log(jj),globals.FloatingPointFormat),',',num2str(fig.YLim(1)+TickMinorLengthY,globals.FloatingPointFormat),')'];
                            end
                        end
                    end
                end
            end
        end
        if strcmp(fig.YMinorTick,'on') && fig.yLog==1 % YMinorTick and log
            if length(fig.YTick)>1, % at least two points
                if fig.YTick(2)-fig.YTick(1) == 1 % one decade per point
                    StrFile{end+1}='% y-Minor Ticks:';
                    for tck=1:length(fig.YTick)-1, % y=const
                        for jj=1:length(subPoints_log),
                            if fig.YTick(tck)+subPoints_log(jj) <= fig.YLim(2),
                                StrFile{end+1}=['\psline[linewidth=\AxesLineWidth,linecolor=GridColor](',num2str(fig.XLim(1),globals.FloatingPointFormat),',',num2str(fig.YTick(tck)+subPoints_log(jj),globals.FloatingPointFormat),')(',num2str(fig.XLim(1)+TickMinorLengthX,globals.FloatingPointFormat),',',num2str(fig.YTick(tck)+subPoints_log(jj),globals.FloatingPointFormat),')'];
                            end
                        end
                    end
                end
            end
        end
        if ~isempty(ahSecond), % additional right side
            subPoints_log2=yAxisCenversionFac*subPoints_log;
            if strcmp(fig.YMinorTick2,'on') && fig.yLog2==1 % YMinorTick and log
                if length(fig.YTick2)>1, % at least two points
                    if fig.YTick2(2)-fig.YTick2(1) == 1 % one decade per point
                        StrFile{end+1}='% y-Minor Ticks Second Axis:';
                        for tck=1:length(fig.YTick2)-1, % y=const
                            for jj=1:length(subPoints_log),
                                StrFile{end+1}=['\psline[linewidth=\AxesLineWidth,linecolor=GridColor](',num2str(fig.XLim(2)-TickMinorLengthX,globals.FloatingPointFormat),',',num2str(fig.YTick2Tex(tck)+subPoints_log2(jj),globals.FloatingPointFormat),')(',num2str(fig.XLim(2),globals.FloatingPointFormat),',',num2str(fig.YTick2Tex(tck)+subPoints_log2(jj),globals.FloatingPointFormat),')'];
                            end
                        end
                    end
                end
            end
        end
        
        
        
    end
    
    % draw labels
    StrFile{end+1}='';
    StrFile{end+1}=['{ ',globals.FontSizeTickLabels,' % FontSizeTickLabels'];
    StrFile{end+1}='% Draw x-Labels: ----';
    for tck=1:length(fig.XTick),
        StrFile{end+1}=['\rput[t](',num2str(fig.XTick(tck),globals.FloatingPointFormat),',',num2str(fig.YLim(1)-TickMajorLengthY,globals.FloatingPointFormat),'){',fig.XTickLabelTex{tck},'}'];
    end
    StrFile{end+1}='% Draw y-Labels: ----';
    for tck=1:length(fig.YTick),
        StrFile{end+1}=['\rput[r](',num2str(fig.XLim(1)-TickMajorLengthX,globals.FloatingPointFormat),',',num2str(fig.YTick(tck),globals.FloatingPointFormat),'){',fig.YTickLabelTex{tck},'}'];
    end
    if ~isempty(ahSecond), % additional right side
        StrFile{end+1}='% Draw y-Labels Second Axis: ----';
        for tck=1:length(fig.YTick2),
            StrFile{end+1}=['\rput[l](',num2str(fig.XLim(2)+TickMajorLengthX,globals.FloatingPointFormat),',',num2str(fig.YTick2Tex(tck),globals.FloatingPointFormat),'){',fig.YTickLabelTex2{tck},'}'];
        end
    end
    StrFile{end+1} = '} % End FontSizeTickLabels'; % end of     FontSizeXYlabel
    

    if strcmp(fig.AxesVisible,'on') && ~globals.AxesBoxEnd, % draw only in this case
        StrFile{end+1}='';
        StrFile{end+1}='% Draw Axes: ----';
        if strcmp(fig.Box,'off'),
            % only left bottom
            StrFile{end+1}=['\psline[linewidth=\AxesLineWidth](',num2str(fig.XLim(1),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')(',num2str(fig.XLim(2),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')'];
            StrFile{end+1}=['\psline[linewidth=\AxesLineWidth](',num2str(fig.XLim(1),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')(',num2str(fig.XLim(1),globals.FloatingPointFormat),',',num2str(fig.YLim(2),globals.FloatingPointFormat),')'];
            if ~isempty(ahSecond), % additional right
                StrFile{end+1}=['\psline[linewidth=\AxesLineWidth](',num2str(fig.XLim(2),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')(',num2str(fig.XLim(2),globals.FloatingPointFormat),',',num2str(fig.YLim(2),globals.FloatingPointFormat),')'];
            end
        else
            tmpstr = '\psframe[linewidth=\AxesLineWidth,dimen=middle](#1#,#2#)(#3#,#4#)';
            StrFile{end+1}=strfill(tmpstr,{fig.XLim(1),fig.YLim(1),fig.XLim(2),fig.YLim(2)});
        end
    end
    
    
    %% LABELS FOR X- AND Y- AXIS AND TITLE
    % XLabel
    wrstr={};
    xlh = get(fig.axh,'XLabel');
    linestr = '\rput[b](#1#,#2#){';
    xlabelX = sum(fig.XLim)/2;
    xlabelY = figy1; % align at bottom
    linestr = strfill(linestr,{xlabelX,xlabelY});
    wrstr = [wrstr,{linestr}];
    wrstr = [wrstr,{'\begin{tabular}{c}'}];
    xlab = get(xlh,'String');
    StrFile{end+1}=''; % add empty line
    StrFile{end+1}=['{ ',globals.FontSizeXYlabel,' % FontSizeXYlabel'];
    if xlab
        StrFile{end+1}='% x-Label: ----';
        for labstr = xlab'
            wrstr = [wrstr,{strcat((labstr(:,1))','\\')}];
        end;
        wrstr = [wrstr,{'\end{tabular}','}'}];
        for str=1:length(wrstr),
            StrFile{end+1}=wrstr{str};
        end
    end;
    
    % YLabel
    wrstr={};
    ylh = get(fig.axh,'YLabel');
    linestr = '\rput[t]{90}(#1#,#2#){';
    ylabelX = figx1; % align left
    ylabelY = sum(fig.YLim)/2;
    linestr = strfill(linestr,{ylabelX,ylabelY});
    wrstr = [wrstr,{linestr}];
    wrstr = [wrstr,{'\begin{tabular}{c}'}];
    ylab = get(ylh,'String');
    if ylab;
        StrFile{end+1}=''; % add empty line
        StrFile{end+1}='% y-Label: ----';
        for labstr = ylab'
            wrstr = [wrstr,{strcat((labstr(:,1))','\\')}];
        end;
        wrstr = [wrstr,{'\end{tabular}'},{'}'}];
        for str=1:length(wrstr),
            StrFile{end+1}=wrstr{str};
        end
    end;
    
    % Second YLabel
    if ~isempty(ahSecond), % additional right side
        wrstr={};
        ylh = get(fig.axh2,'YLabel');
        linestr = '\rput[b]{90}(#1#,#2#){';
        ylabelX = figx2; % align right
        ylabelY = sum(fig.YLim)/2;
        linestr = strfill(linestr,{ylabelX,ylabelY});
        wrstr = [wrstr,{linestr}];
        wrstr = [wrstr,{'\begin{tabular}{c}'}];
        ylab = get(ylh,'String');
        if ylab;
            StrFile{end+1}=''; % add empty line
            StrFile{end+1}='% y-Label Second Axis: ----';
            for labstr = ylab'
                wrstr = [wrstr,{strcat((labstr(:,1))','\\')}];
            end;
            wrstr = [wrstr,{'\end{tabular}'},{'}'}];
            for str=1:length(wrstr),
                StrFile{end+1}=wrstr{str};
            end
        end
    end
    
    % Title
    th = get(fig.axh,'Title');
    wrstr={};
    linestr = '\rput[t](#1#,#2#){';
    titleX = sum(fig.XLim)/2;
    titleY = figy2; % align at top
    linestr = strfill(linestr,{titleX,titleY});
    wrstr = [wrstr,{linestr}];
    wrstr = [wrstr,{'\begin{tabular}{c}'}];
    titlab = get(th,'String');
    if titlab
        StrFile{end+1}=''; % add empty line
        StrFile{end+1}='% Title: ----';
        for labstr = titlab'
            wrstr = [wrstr,{strcat((labstr(:,1))','\\')}];
        end;
        wrstr = [wrstr,{'\end{tabular}'},{'}'}];
        for str=1:length(wrstr),
            StrFile{end+1}=wrstr{str};
        end
    end;
    
    StrFile{end+1} = '} % End FontSizeXYlabel'; % end of     FontSizeXYlabel

end % axes visible


if ischar(globals.AxesLineWidth)
    StrFile={['\providelength{\AxesLineWidth}       \setlength{\AxesLineWidth}{',globals.AxesLineWidth,'}%'],StrFile{1:end}};
else
    StrFile={['\providelength{\AxesLineWidth}       \setlength{\AxesLineWidth}{',num2str(globals.AxesLineWidth),'pt}%'],StrFile{1:end}};
end
StrFile={'% Global Parameters that can be changed:',StrFile{1:end}};

%% WRITE DATA TO FILE
lshandles = get(fig.axh,'Children');
lshandles = double(find(handle(lshandles),'-class','graph2d.lineseries'));
lshandles = flipud(lshandles); % flip lshandles in order to get the correct plot order
kk = 0;
for lsh=lshandles'; % scan all data lines
    kk = kk + 1;
    
    xd = get(lsh,'XData');
    yd = get(lsh,'YData');
    
    if length(xd)==1, % a single point, double for visualization
        xd=[xd, xd]; yd=[yd, yd];
    end
    
    % cut data that are out of bounds and append them as extra lines:
    [Pcutt, Psnipp] = cut_data (xd, yd, fig.XLim, fig.YLim);

   
    xd_cut = {};
    yd_cut = {};
    mark_on = [];
    
    % append snippets
    for snp = 1 : length(Psnipp)
        xd_cut{snp} = Psnipp{snp}(:,1);
        yd_cut{snp} = Psnipp{snp}(:,2);
        mark_on(snp) = 0;
    end
    
    % append valid lines
    for ln = 1 : length(Pcutt)        
        if size(Pcutt{ln}, 1) == 1,  % a single point, double for visualization
            xd_cut{end+1} = [Pcutt{ln}(:,1); Pcutt{ln}(:,1)];
            yd_cut{end+1} = [Pcutt{ln}(:,2); Pcutt{ln}(:,2)];
        else
            xd_cut{end+1} = Pcutt{ln}(:,1);
            yd_cut{end+1} = Pcutt{ln}(:,2);            
        end
        mark_on(end+1) = 1;
    end

    % copy back:
    xd = xd_cut;
    yd = yd_cut;
     
    % get the markers
    markerStr=convert_Marker(get(lsh,'Marker'));
    markerStr=markerStr{1}; % no cell array
    

    % get the linestyle
    linestyle = convert_lineStyle(get(lsh,'LineStyle'));
    linestyle = linestyle{1}; % no cell array
    

    % check for patch
    if isempty(fig.FaceColor{kk}), % line, no patch
        isPatch = false;
        StrFile{end+1}='';
        cname   = ['color',num2str(lsh)];
        cnameFC = [];
        colstr  = ['\newrgbcolor{',cname,'}{',num2str(get (lsh,'Color')),'}'];
        StrFile{end+1}='% New Line DATA: ----';
        StrFile{end+1}=colstr;
        
    else
        isPatch = true;
        
        % put all snippets together in order to obtain a single cut patch again
        if ~isempty(xd),
            xdd = xd{1};
            ydd = yd{1};
            for jj=2:length(xd)-1, % all cut data are written
                xdd = [xdd; xd{jj}];
                ydd = [ydd; yd{jj}]; % save current cut data
            end
            xd = {xdd};
            yd = {ydd};

            StrFile{end+1}='';
            Fcol = fig.FaceColor{kk};
            Lcol = fig.LineColor{kk};
            cname   = ['Lcolor',num2str(lsh)];
            cnameFC = ['Fcolor',num2str(lsh)];
            StrFile{end+1}='% New Patch DATA: ----';
            colstr  = ['\newrgbcolor{',cname,'}{',num2str(fig.LineColor{kk}),'}'];
            StrFile{end+1}=colstr;
            colstr  = ['\newrgbcolor{',cnameFC,'}{',num2str(fig.FaceColor{kk}),'}'];
            StrFile{end+1}=colstr;

        end

    end
    
    % write data
    for jj=1:length(xd), % all cut data are written
        xdd=xd{jj}'; ydd=yd{jj}'; % save current cut data
       
        if length (xdd) > 1, % at least two points
            
            % reduce number of points for the current line
            if strcmp(markerStr,'none'), % do NOT reduce for lines with markers
                xydata=reduce_line([xdd;ydd],globals.PointReduceFactor);
                xdd=xydata(1,:);
                ydd=xydata(2,:);
            end
            
            if isPatch % in this case patches are plotted
                
                StrFile{end+1}='\savedata{\mydata}[{';
                tmpstr=''; % init
                for str=1:length(xdd),
                    tmpstr = [tmpstr,'{',num2str(xdd(str),globals.FloatingPointFormat),',',num2str(ydd(str),globals.FloatingPointFormat),'},'];
                    if mod(str,5)==0, % new line all 5 data pair
                        StrFile{end+1}=tmpstr;
                        tmpstr='';
                    end
                end
                if not(isempty(tmpstr)), % last line
                    tmpstr(end)=''; % delete last ','
                    StrFile{end+1}=tmpstr;
                end
                StrFile{end+1}='}]';
                                
                tmpstr = ['\dataplot[fillstyle=solid,fillcolor=#2#,plotstyle=line,linestyle=#1#,linewidth=\EdgeLineWidth,linecolor=#3#]{\mydata}'];
                wrstr = strfill(tmpstr,{linestyle,cnameFC,cname});
                StrFile{end+1}=wrstr;
                
            else % line
                if strcmp(markerStr,'none'), % do not show markers
                    tmpstr = ['\psline[plotstyle=line,linejoin=1,linestyle=#1#,linewidth=\LineWidth,linecolor=#2#]'];
                else                         % show markers
                    if mark_on(jj),
                        showpoints = 'true';
                    else
                        showpoints = 'false';
                    end
                    tmpstr = ['\psline[plotstyle=line,linejoin=1,showpoints=',showpoints,',dotstyle=',markerStr,',dotsize=\MarkerSize,linestyle=#1#,linewidth=\LineWidth,linecolor=#2#]'];
                end
                wrstr = strfill(tmpstr,{linestyle,cname});
                StrFile{end+1}=wrstr;
                
                tmpstr='';
                for str=1:length(xdd),
                    tmpstr = [tmpstr,'(',num2str(xdd(str),globals.FloatingPointFormat),',',num2str(ydd(str),globals.FloatingPointFormat),')'];
                    if mod(str,5)==0, % new line all 5 data pair
                        StrFile{end+1}=tmpstr;
                        tmpstr='';
                    end
                end
                if not(isempty(tmpstr)), % last line
                    StrFile{end+1}=tmpstr;
                end
                
                
            end
            
        end
    end
end

%% Axes Box
if strcmp(fig.AxesVisible,'on') && globals.AxesBoxEnd, % draw only in this case
    StrFile{end+1}='';
    StrFile{end+1}='% Draw Axes: ----';
    if strcmp(fig.Box,'off'),
        % only left bottom
        StrFile{end+1}=['\psline[linewidth=\AxesLineWidth](',num2str(fig.XLim(1),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')(',num2str(fig.XLim(2),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')'];
        StrFile{end+1}=['\psline[linewidth=\AxesLineWidth](',num2str(fig.XLim(1),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')(',num2str(fig.XLim(1),globals.FloatingPointFormat),',',num2str(fig.YLim(2),globals.FloatingPointFormat),')'];
        if ~isempty(ahSecond), % additional right
            StrFile{end+1}=['\psline[linewidth=\AxesLineWidth](',num2str(fig.XLim(2),globals.FloatingPointFormat),',',num2str(fig.YLim(1),globals.FloatingPointFormat),')(',num2str(fig.XLim(2),globals.FloatingPointFormat),',',num2str(fig.YLim(2),globals.FloatingPointFormat),')'];
        end
    else
        tmpstr = '\psframe[linewidth=\AxesLineWidth,dimen=middle](#1#,#2#)(#3#,#4#)';
        StrFile{end+1}=strfill(tmpstr,{fig.XLim(1),fig.YLim(1),fig.XLim(2),fig.YLim(2)});
    end
end



%% LEGEND
[legh,~,plot_h,text_strings] = legend(fig.axh); % get legend information
if ~isempty(legh)
    ts=globals.LegendBoxTeXcode; % first line for the legend box
    % ts={'\psshadowbox[framesep=0pt]{\psframebox*{\begin{tabular}{l}'};
    rowstr='\Rnode{a#1#}{\hspace*{0.0ex}} \hspace*{0.7cm} \Rnode{a#2#}{~~#3#} \\';
    llstr={};
    if numPatches > length(plot_h),
        kstart = 1;
        nPatchs = 0;
    else
        kstart  = numPatches+1;
        nPatchs = numPatches;
    end
    for k=kstart:length(plot_h)% do not account for the patches at the beginning
        markerStr=convert_Marker(get(plot_h(k),'Marker'));
        markerStr=markerStr{1}; % no cell array
        if strcmp(markerStr,'none'), % do not show markers
            ncstr='\ncline[linestyle=#1#,linewidth=\LineWidth,linecolor=#4#]{a#2#}{a#3#}';
        else
            ncstr=['\ncline[linestyle=#1#,linewidth=\LineWidth,linecolor=#4#]{a#2#}{a#3#} \ncput{\psdot[dotstyle=',markerStr,',dotsize=\MarkerSize,linecolor=#4#]}'];
        end
        tabstr = strfill(rowstr,{num2str(2*k-1),num2str(2*k),text_strings(k-nPatchs)});
        ts=[ts,tabstr];
        linestyle = convert_lineStyle(get(plot_h(k),'LineStyle'));
        linestyle = linestyle{1}; % convert from cell to chars
        cname=['color',num2str(plot_h(k))];
        leglinestr = strfill(ncstr,{linestyle,num2str(2*k-1),num2str(2*k),cname});
        llstr = [llstr,{leglinestr}];
    end;
    ts=[ts,{'\end{tabular}}'}];
    
    % legend's position
    set(legend_org,'Units','normalized')
    lp = get(legend_org,'Position'); % position normalized
    % convert relative figure positions to x-and y-values
    set(ah,'Units','normalized');
    pos=get(ah,'Position');
    XY1=pos(1:2); Width=pos(3:4);
    legendposX=((lp(1)+lp(3)/2)-XY1(1))./Width(1)*(fig.XLim(2)-fig.XLim(1))+fig.XLim(1); % normalized axes units
    legendposY=((lp(2)+lp(4)/2)-XY1(2))./Width(2)*(fig.YLim(2)-fig.YLim(1))+fig.YLim(1);
    
    switch get(legh,'Location'),
        case 'NorthWest'
            tmpstr = '\rput[tl](#1#,#2#){%';
            legendposX = fig.XLim(1)+TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = fig.YLim(2)-TickMajorLengthY*globals.LegendOffsetFac;
        case 'SouthWest'
            tmpstr = '\rput[bl](#1#,#2#){%';
            legendposX = fig.XLim(1)+TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = fig.YLim(1)+TickMajorLengthY*globals.LegendOffsetFac;
        case 'SouthEast'
            tmpstr = '\rput[br](#1#,#2#){%';
            legendposX = fig.XLim(2)-TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = fig.YLim(1)+TickMajorLengthY*globals.LegendOffsetFac;
        case 'NorthEast'
            tmpstr = '\rput[tr](#1#,#2#){%';
            legendposX = fig.XLim(2)-TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = fig.YLim(2)-TickMajorLengthY*globals.LegendOffsetFac;
        case 'West'
            tmpstr = '\rput[l](#1#,#2#){%';
            legendposX = fig.XLim(1)+TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = (fig.YLim(1) + fig.YLim(2)) / 2;
        case 'South'
            tmpstr = '\rput[b](#1#,#2#){%';
            legendposX = (fig.XLim(1) + fig.XLim(2)) / 2;
            legendposY = fig.YLim(1)+TickMajorLengthY*globals.LegendOffsetFac;
        case 'East'
            tmpstr = '\rput[r](#1#,#2#){%';
            legendposX = fig.XLim(2)-TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = (fig.YLim(1) + fig.YLim(2)) / 2;
        case 'North'
            tmpstr = '\rput[t](#1#,#2#){%';
            legendposX = (fig.XLim(1) + fig.XLim(2)) / 2;
            legendposY = fig.YLim(2)-TickMajorLengthY*globals.LegendOffsetFac;
        case 'NorthEastOutside'
            tmpstr = '\rput[lt](#1#,#2#){%';
            legendposX = fig.XLim(2)+TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = fig.YLim(2);
        case 'EastOutside'
            tmpstr = '\rput[l](#1#,#2#){%';
            legendposX = fig.XLim(2)+TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = (fig.YLim(1) + fig.YLim(2)) / 2;
        case 'SouthEastOutside'
            tmpstr = '\rput[lb](#1#,#2#){%';
            legendposX = fig.XLim(2)+TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = fig.YLim(1);
        case 'NorthWestOutside'
            tmpstr = '\rput[rt](#1#,#2#){%';
            legendposX = figx1-TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = fig.YLim(2);
        case 'WestOutside'
            tmpstr = '\rput[r](#1#,#2#){%';
            legendposX = figx1-TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = (fig.YLim(1) + fig.YLim(2)) / 2;
        case 'SouthWestOutside'
            tmpstr = '\rput[rb](#1#,#2#){%';
            legendposX = figx1-TickMajorLengthX*globals.LegendOffsetFac;
            legendposY = fig.YLim(1);
        otherwise
            tmpstr = '\rput(#1#,#2#){%';
    end
    wrstr = strfill(tmpstr,{legendposX,legendposY});
    StrFile{end+1}=''; StrFile{end+1}='% Legend: ----';
    StrFile{end+1}=['{ ',globals.FontSizeLegend,' % FontSizeLegend'];
    StrFile{end+1}=wrstr;
    
    tmpstr2 = [ts,llstr,{'}%'},{'}%'}];
    for str=1:length(tmpstr2)
        StrFile{end+1}=tmpstr2{str};
    end;
    StrFile{end+1} = '} % End FontSizeLegend'; % end of     FontSizeLegend
    
end;

%% Colorbar
if ~isempty(ColBar),
    
    if strcmp(ColBar.Location, 'EastOutside'), % only location supported right now
        
        StrFile{end+1}=''; StrFile{end+1}='% Colorbar: ----';
        
        % Position of the colorbar
        
        lp = ColBar.Position; % position normalized
        % convert relative figure positions to x-and y-values
        set(ah,'Units','normalized');
        pos=get(ah,'Position');
        XY1=pos(1:2); Width=pos(3:4);
        ColBarPosYB=fig.YLim(1);
        ColBarPosYT=fig.YLim(2);
        if isempty(globals.ColorBarWidth_cm)
            
            % (FROM MATLAB)
            ColBarPosXL=((lp(1))      -XY1(1))./Width(1)*(fig.XLim(2)-fig.XLim(1))+fig.XLim(1); % normalized axes units
            ColBarPosXR=((lp(1)+lp(3))-XY1(1))./Width(1)*(fig.XLim(2)-fig.XLim(1))+fig.XLim(1);
            % ColBarPosYB=((lp(2))      -XY1(2))./Width(2)*(fig.YLim(2)-fig.YLim(1))+fig.YLim(1);
            % ColBarPosYT=((lp(2)+lp(4))-XY1(2))./Width(2)*(fig.YLim(2)-fig.YLim(1))+fig.YLim(1);
            
        elseif length(globals.ColorBarWidth_cm) == 1, % scalar width
        
            ColBarPosXL=((lp(1))-XY1(1))./Width(1)*(fig.XLim(2)-fig.XLim(1))+fig.XLim(1); % normalized axes units
            ColBarPosXR = ColBarPosXL+globals.ColorBarWidth_cm/globals.plotwidth*DeltaXLim;
            
        elseif length(globals.ColorBarWidth_cm) == 2, % [width, offset]
            
            ColBarPosXL = fig.XLim(2) + globals.ColorBarWidth_cm(2)/globals.plotwidth*DeltaXLim; % normalized axes units
            ColBarPosXR = ColBarPosXL + globals.ColorBarWidth_cm(1)/globals.plotwidth*DeltaXLim;
           
        else
            error('Wrong input format of globals.ColorBarWidth_cm!')
        end
        
        cmap = ColBar.cmap;
        
        % Draw all color patches
        numColPatches = ColBar.CLim(2) - ColBar.CLim(1);
        % numColPatches = size(cmap, 1) % default is 64
        dY = (ColBarPosYT - ColBarPosYB) / numColPatches;
        for pt = 1 : numColPatches,
            randName = num2str(rand,8); % random name
            cname = ['colB',num2str(pt),'_',randName(3:6)];
            StrFile{end+1} = ['\newrgbcolor{',cname,'}{', ...
                num2str(cmap(pt,1)),' ', ...
                num2str(cmap(pt,2)),' ', ...
                num2str(cmap(pt,3)),'}'];
            tmpstr = '\psframe[linewidth=0.1pt,linestyle=solid,linecolor=#1#,fillstyle=solid,fillcolor=#1#,dimen=middle](#2#,#3#)(#4#,#5#)';
            StrFile{end+1}=strfill(tmpstr,{cname, ...
                ColBarPosXL,ColBarPosYB+(pt-1)*dY,ColBarPosXR,ColBarPosYB+pt*dY});
        end
        % Draw border
        tmpstr = '\pspolygon[linewidth=\AxesLineWidth](#1#,#2#)(#3#,#2#)(#3#,#4#)(#1#,#4#)(#1#,#2#)';
        StrFile{end+1}=strfill(tmpstr,{ColBarPosXL,ColBarPosYB,ColBarPosXR,ColBarPosYT});
        
        
        % Convert positions of ticks
        numTicks = length(ColBar.YTick);
        for tck = 1 : numTicks,
            Yact(tck) = ColBarPosYB + (ColBar.YTick(tck) - ColBar.YLim(1)) * ...
                (ColBarPosYT - ColBarPosYB) / (ColBar.YLim(2) - ColBar.YLim(1));
        end
        
        % draw ticks
        TickColor = ColBar.YColor;
        randName = num2str(rand,8); % random name
        TickColorName = ['ColBarTickCol_',randName(3:6)];
        StrFile{end+1}='';
        StrFile{end+1}='% Draw Colorbar y-Ticks: ----';
        tmpstr = '\newrgbcolor{#1#}{#2# #3# #4#}';
        StrFile{end+1}=strfill(tmpstr,{ ...
            TickColorName, TickColor(1), TickColor(2), TickColor(3)});
        ColBarWidth = ColBarPosXR - ColBarPosXL;
        for tck=1:numTicks, % y=const
            % left tick:
            tmpstr = '\psline[linewidth=\AxesLineWidth,linecolor=#4#](#1#,#2#)(#3#,#2#)';
            StrFile{end+1}=strfill(tmpstr, ...
                {ColBarPosXL,Yact(tck),ColBarPosXL+ColBarWidth/6, ...
                TickColorName});
            % right tick:
            tmpstr = '\psline[linewidth=\AxesLineWidth,linecolor=#4#](#1#,#2#)(#3#,#2#)';
            StrFile{end+1}=strfill(tmpstr, ...
                {ColBarPosXR,Yact(tck),ColBarPosXR-ColBarWidth/6, ...
                TickColorName});
        end
        
        % draw labels
        StrFile{end+1}='';
        StrFile{end+1}=['{ ',globals.FontSizeTickLabels,' % FontSizeTickLabels'];
        StrFile{end+1}='% Draw Colorbar y-Labels: ----';
        labels_out = convert_TickLabels2Cell (ColBar.YTick, ColBar.YTickLabel);
        for tck=1:numTicks,
            tmpstr = '\rput[l](#1#,#2#){$#3#$}';
            StrFile{end+1}=strfill(tmpstr, ...
                {ColBarPosXR+TickMajorLengthX,Yact(tck),labels_out{tck}});
        end
        StrFile{end+1} = '} % End FontSizeTickLabels'; % end of     FontSizeXYlabel
        
        
        % drax y-label
        StrFile{end+1}=''; % add empty line
        StrFile{end+1}=['{ ',globals.FontSizeXYlabel,' % FontSizeXYlabel'];
        wrstr={};
        linestr = '\rput[b]{90}(#1#,#2#){';
        % align right and use the same offset the y label on the left hand
        % side:
        % ColBarPosXR = ColBarPosXL + globals.ColorBarWidth_cm(1)/globals.plotwidth*DeltaXLim;
        
        % test if right tight inset is larger than left tight inset +
        % colorbar width
        if globals.tightInset_cm(3) > ...
                (globals.tightInset_cm(1) + sum(globals.ColorBarWidth_cm)),
            ylabelX = figx2; % right
        else
            ylabelX = ColBarPosXR + globals.tightInset_cm(1)/globals.plotwidth*DeltaXLim;
        end
        ylabelY = sum(fig.YLim)/2;
        linestr = strfill(linestr,{ylabelX,ylabelY});
        wrstr = [wrstr,{linestr}];
        wrstr = [wrstr,{'\begin{tabular}{c}'}];
        ylab = get(ColBar.YLabel,'String');
        if ylab;
            StrFile{end+1}='% y-Label Colorbar: ----';
            for labstr = ylab'
                wrstr = [wrstr,{strcat((labstr(:,1))','\\')}];
            end;
            wrstr = [wrstr,{'\end{tabular}'},{'}'}];
            for str=1:length(wrstr),
                StrFile{end+1}=wrstr{str};
            end
        end
        StrFile{end+1} = '} % End FontSizeXYlabel'; % end of  FontSizeXYlabel
        
        
        
    else
        
        warning('The only supported location of the Colorbar is EastOutside!')
        
    end
    
end





%% Textarrows
if ~isempty(anhar),
    
    if ~ischar (globals.TextArrowLineWidth)
        TextArrowLineWidth = [num2str(globals.TextArrowLineWidth),'pt'];
    else
        TextArrowLineWidth = globals.TextArrowLineWidth;
    end
    
    StrFile{end+1}='';
    StrFile{end+1}=['{ ',globals.FontSizeAnnotations,' % FontSizeAnnotations'];
    StrFile{end+1}='% Text Arrows: ----';
    for str=1:length(fig.TextarrowsStr),
        wrstr ={}; % init
        wrstr = [wrstr,['% Number ',num2str(str),':']];
        cname=['color',num2str(anhar(str))]; % name is random, just num2str of the handle
        colstr=['\newrgbcolor{',cname,'}{',num2str(fig.TextarrowsColor{str}),'}'];
        wrstr = [wrstr,colstr];
        annostr = fig.TextarrowsStr{str};        
        annoX = fig.TextarrowsX(str,:);
        annoY = fig.TextarrowsY(str,:);
        lstTA=convert_lineStyle(fig.TextarrowsLineStyle{str});
        linestr = ['\psline[linestyle=',lstTA{1},',linewidth=', ...
            TextArrowLineWidth,',linecolor=#5#,arrowsize=1.5pt 3,arrowlength=2,arrowinset=0.3]{->}(#1#,#2#)(#3#,#4#)'];
        linestr = strfill(linestr,{annoX(1),annoY(1),annoX(2),annoY(2),cname});
        wrstr = [wrstr,linestr];
        if ~isempty(annostr), % do this only if a string is present
            linestr = '\uput{0pt}[#1#](#2#,#3#){%';
            linestr = strfill(linestr,{fig.TextarrowsRot(str),annoX(1),annoY(1)});
            wrstr = [wrstr,linestr];
            if length(annostr)==1, % text is just one line, don't use a tabular
                if globals.TextArrowTextColor==1, % adapt the color
                    linestr='\psframebox*[framesep=1pt]{\textcolor{#2#}{#1#}}}';
                    linestr = strfill(linestr,{annostr,cname});
                else
                    linestr='\psframebox*[framesep=1pt]{#1#}}';
                    linestr = strfill(linestr,{annostr});
                end
                wrstr = [wrstr,linestr];
            else
                if globals.TextArrowTextColor==1, % adapt the color
                    if strcmp(fig.TextarrowsHorizontalAlignment, 'left'),
                        al = 'l';
                    elseif strcmp(fig.TextarrowsHorizontalAlignment, 'right'),
                        al = 'r';
                    else
                        al = 'c';
                    end
                    linestr =['\psframebox*[framesep=1pt]{\textcolor{#1#}{\begin{tabular}{@{}',al,'@{}}'];
                    linestr = strfill(linestr,{cname});
                    endstr  =  ['\end{tabular}','}}}'];
                else
                    linestr='\psframebox*[framesep=1pt]{\begin{tabular}{@{}c@{}}';
                    endstr =  ['\end{tabular}','}}'];
                end
                wrstr = [wrstr,linestr];
                if ~iscell (annostr), annostr = {annostr}; end
                for jj=1:length(annostr),
                    wrstr = [wrstr,strcat(annostr{jj},'\\[-0.3ex]')]; % [-0.3ex] negative space
                end
                wrstr = [wrstr,endstr];
            end
        end
        for ll=1:length(wrstr), % write everything to file
            StrFile{end+1}=wrstr{ll};
        end
    end
    StrFile{end+1} = '} % End FontSizeAnnotations'; % end of     FontSizeLegend
end


%% Text boxes
if ~isempty(anhtb),
    
    if ~ischar (globals.TextBoxLineWidth)
        TextBoxLineWidth = [num2str(globals.TextBoxLineWidth),'pt'];
    else
        TextBoxLineWidth = globals.TextBoxLineWidth;
    end
    
    StrFile{end+1}='';
    StrFile{end+1}=['{ ',globals.FontSizeAnnotations,' % FontSizeAnnotations'];
    StrFile{end+1}='% Text Boxes: ----';
    for str=1:length(fig.TextBoxesStr),
        
        % convert current linestyle of the box
        lstyle = convert_lineStyle(fig.TextBoxesLineStyle{str});

        wrstr ={}; % init
        wrstr = [wrstr,['% Number ',num2str(str),':']];
        cname=['color',num2str(anhtb(str))]; % name is random, just num2str of the handle
        colstr=['\newrgbcolor{',cname,'}{',num2str(fig.TextBoxesColor{str}),'}'];
        wrstr = [wrstr,colstr];
        annostr = fig.TextBoxesStr{str};
        annoX = fig.TextBoxesX(str,:);
        annoY = fig.TextBoxesY(str,:);
        linestr = '\uput{0pt}[0](#1#,#2#){%';
        linestr = strfill(linestr,{annoX(1),annoY(1)});
        wrstr = [wrstr,linestr];
         if length(annostr)==1, % text is just one line, don't use a tabular
            if globals.TextBoxTextColor==1, % adapt the color
                linestr='\psframebox[framesep=1pt,fillstyle=solid,linestyle=#3#,linewidth=#4#]{\textcolor{#2#}{#1#}}}';
                linestr = strfill(linestr,{annostr,cname,lstyle,TextBoxLineWidth});
            else
                linestr='\psframebox[framesep=1pt,fillstyle=solid,linestyle=#2#,linewidth=#3#]{#1#}}';
                linestr = strfill(linestr,{annostr,lstyle,TextBoxLineWidth});
            end
            wrstr = [wrstr,linestr];
        else
            if globals.TextBoxTextColor==1, % adapt the color
                if strcmp(fig.TextBoxesHorizontalAlignment, 'left'),
                    al = 'l';
                elseif strcmp(fig.TextBoxesHorizontalAlignment, 'right'),
                    al = 'r';
                else
                    al = 'c';
                end
                linestr =['\psframebox[framesep=1pt,fillstyle=solid,linestyle=#2#,linewidth=#3#]{\textcolor{#1#}{\begin{tabular}{@{}',al,'@{}}'];
                linestr = strfill(linestr,{cname,lstyle,TextBoxLineWidth});
                endstr  =  ['\end{tabular}','}}}'];
            else
                linestr='\psframebox[framesep=1pt,fillstyle=solid,linestyle=#1#,linewidth=#2#]{\begin{tabular}{@{}c@{}}';
                linestr = strfill(linestr,{lstyle,TextBoxLineWidth});
                endstr =  ['\end{tabular}','}}'];
            end
            wrstr = [wrstr,linestr];
            if ~iscell (annostr), annostr = {annostr}; end
            for jj=1:length(annostr),
                wrstr = [wrstr,strcat(annostr{jj},'\\[-0.3ex]')]; % [-0.3ex] negative space
            end
            wrstr = [wrstr,endstr];
        end
        for ll=1:length(wrstr), % write everything to file
            StrFile{end+1}=wrstr{ll};
        end
    end
    StrFile{end+1} = '} % End FontSizeAnnotations'; % end of     FontSizeLegend
end


%% Convergence Triangle

% globals.TriangleOrder = []; % vector with orders that are shown, if empty, no convergence triangle
% globals.TriangleColor = [0.4 0.4 0.4]; % rgb vector
% globals.TriangleLxRel = 0.3;

if ~isempty(globals.TriangleOrder)
    ConTriposX = fig.XLim(1)+TickMajorLengthX*globals.LegendOffsetFac;
    ConTriposY = fig.YLim(1)+TickMajorLengthY*globals.LegendOffsetFac;
    bordersX   = fig.XLim;
    StartPoint = [ConTriposX, ConTriposY];
    numOrders  = length(globals.TriangleOrder);
    LineStyle  = 'solid';
    labelPos   = 0.2;
    
    StrFile{end+1}='';
    StrFile{end+1}=['{ ',globals.FontSizeAnnotations,' % FontSizeAnnotations'];
    StrFile{end+1}='% Convergence Triangle: ----';
    StrFile{end+1}=['\newrgbcolor{TriangleColor}{',num2str(globals.TriangleColor),'}'];
    
    % Nodes
    LB=StartPoint;
    RB=StartPoint+globals.TriangleLxRel*(bordersX(2)-bordersX(1))*[1 0];
    StrFile{end+1}=['\pnode(',num2str(LB(1)),',',num2str(LB(2)),'){LB}'];
    StrFile{end+1}=['\pnode(',num2str(RB(1)),',',num2str(RB(2)),'){RB}'];
    for or = 1 : numOrders,
        LT{or} = StartPoint + globals.TriangleLxRel * ...
            (bordersX(2)-bordersX(1)) * globals.TriangleOrder(or) * [0 1];
        StrFile{end+1}=['\pnode(',num2str(LT{or}(1)),',',num2str(LT{or}(2)),'){LT',num2str(or),'}'];
    end
    
    % Lines
    StrFile{end+1}=['\ncline[linecolor=TriangleColor,linestyle=',LineStyle,',linewidth=\AxesLineWidth,arrows=c-c]{LB}{RB}'];
    StrFile{end+1}=['\ncline[linecolor=TriangleColor,linestyle=',LineStyle,',linewidth=\AxesLineWidth,arrows=c-c]{LB}{LT',num2str(numOrders),'}'];
    for or = 1 : numOrders,
        if ~isempty(globals.TriangleStr), % extra string
            str = globals.TriangleStr{or};
        else
            str = sprintf('$%.2f$',globals.TriangleOrder(or));
        end        
        StrFile{end+1}=['\ncline[linecolor=TriangleColor,linestyle=',LineStyle,',linewidth=\AxesLineWidth,arrows=c-c]{LT',num2str(or),'}{RB}'];
        StrFile{end+1}=['\ncput[npos=',num2str(labelPos),']{\psframebox*[framesep=1pt]{',str,'}}'];
    end
    
    StrFile{end+1} = '} % End FontSizeAnnotations'; % end of     FontSizeLegend
    
end



%%

for cd = 1 : length(globals.AddCode) % append additional external code
    StrFile{end+1} = globals.AddCode{cd};
end

StrFile{end+1}='';
StrFile{end+1}='\end{pspicture}%';

% Add beginning strings
if globals.StandAloneFile,
    
    tempstr={['% Date:  ',datestr(now)],...
        '% This file was created by fig2texps.','',...
        '\documentclass{standalone}', ...
        '\usepackage{ifpdf}', ...
        '\ifpdf% Include pdf wrapper if running pdflatex ',...
        '  \usepackage[pdf]{pstricks}%', ...
        '  \usepackage{auto-pst-pdf}%', ...
        '\else%', ...
        '  \usepackage{pstricks}%', ...
        '\fi%', ...
        '\usepackage{pst-node, pst-plot, pst-circ}', ...
        '\usepackage{moredefs}', ...
        '\begin{document}',''};
        %'\usepackage[%', ...        
        %['    paperwidth = ',num2str(TotalWidth_cm*1.005,8),'cm, % add 0.5%'], ...
        %['    paperheight= ',num2str(TotalHeight_cm*1.005,8),'cm, % add 0.5%'], ...
        %'    lmargin    = -5.25mm, % do not know why this offset is needed.', ...
        %'    rmargin    = 0mm, %', ...
        %'    tmargin    = 0mm, %', ...
        %'    bmargin    = 0mm, %', ...
        %'    ]{geometry}', ...
    StrFile={tempstr{:},StrFile{1:end}};
    StrFile={StrFile{1:end}, '', '\end{document}'};
    
else
    
    tempstr={['% Date:  ',datestr(now)],...
        '% This file was created by fig2texps. Note, that the packages',...
        '% pstricks, pst-node, pst-plot, pst-circ, and moredefs are required.', ...        
        '% Also, ifpdf and auto-pst-pdf packages are require to run pdflatex.', ...
        '% If so, remember to run pdflatex --shell-escape file.tex.', ...
        '% If you want to improve the speed of compilation, consider using standalone package.', ...
        '% A minimal example code could be:','%',...
        '% \documentclass{article}', ...
        '% \usepackage{ifpdf}', ...
        '% \ifpdf% Include pdf wrapper if running pdflatex ',...
        '%   \usepackage[pdf]{pstricks}%', ...
        '%   \usepackage{auto-pst-pdf}%', ...
        '% \else%', ...
        '%  \usepackage{pstricks}%', ...
        '% \fi%', ...
        '% \usepackage{pst-node, pst-plot, pst-circ}', ...
        '% \usepackage{moredefs}', ...
        '% \begin{document}', ...
        '% \input{fig1.tex}', ...
        '% \end{document}','%'};
    StrFile={tempstr{:},StrFile{1:end}};
end

warning(s); % turn on specified warnings again
% PRINT TO FILE

% Check if file already exists
if exist(fullfile(globals.pathname,globals.filename)),
    % Construct a questdlg with three options
    choice = questdlg('File already exists. Replace it?', ...
        'Overwrite file', ...
        'Yes','No','Cancel','No');
    switch choice
        case 'Yes'
            fid = fopen(fullfile(globals.pathname,globals.filename),'wt');
            for str=1:length(StrFile)-1,
                wrstr=StrFile{str};
                fprintf(fid,'%s\n',wrstr);
            end
            fprintf(fid,'%s',StrFile{end}); % last line without 'return'
            fclose(fid);
            disp('File successfully written.')
    end
    
else
    
    fid = fopen(fullfile(globals.pathname,globals.filename),'wt');
    for str=1:length(StrFile)-1,
        wrstr=StrFile{str};
        fprintf(fid,'%s\n',wrstr);
    end
    fprintf(fid,'%s',StrFile{end}); % last line withour 'return'
    fclose(fid);
    disp('File successfully written.')
    
end





% close figure
close(fig.figh)


% -------------------------------------------------------------------------
% END OF MAIN FUNCTION ----------------------------------------------------
% -------------------------------------------------------------------------




%% Additional Functions

% Parse input 
    function parse(varargin)
        i = 1;
        while (i <= length(varargin))
            if (i+1 > length(varargin)), error('Input parameter missed for %s',varargin{i}); end
            if isfield(globals,varargin{i})
                if(strcmp(varargin{i},'filename'))
                    [pathstr, name, ~] = fileparts(fname);
                    if ~isempty(pathstr),
                        globals.pathname = pathstr;
                    end
                    globals.filename = strcat(strrep(name,'.tex',''),'.tex'); % append .tex if not already there
                else
                    globals.(varargin{i}) = varargin{i+1};
                end
            else
                error('Argument %s does not exist',varargin{i});
            end
            i = i+2;
        end
    end

% SUBFUNCTIONS

    function linestyle = convert_lineStyle(lsty)
        % Convert the line style to tex names,
        % lsty is a cell array
        if not(iscell(lsty)), lsty={lsty}; end
        for ii_lst=1:length(lsty)
            if isempty(lsty{ii_lst}),linestyle{ii_lst}='none';
            else
                switch lsty{ii_lst}
                    case 'none'
                        linestyle{ii_lst}='none';
                    case '-'
                        linestyle{ii_lst}='solid';
                    case ':'
                        linestyle{ii_lst}='dashed,dash=2pt 3pt';
                    case '--'
                        linestyle{ii_lst}='dashed';
                    case '-.'
                        linestyle{ii_lst}='dashed,dash=3pt 2pt 1pt 2pt';
                    otherwise
                        linestyle{ii_lst}='solid';
                end
            end
        end
    end


    function marker = convert_Marker(mrk)
        % Convert the line style to tex names,
        % lsty is a cell array
        if not(iscell(mrk)), mrk={mrk}; end
        for ii=1:length(mrk)
            switch mrk{ii}
                case '+',           marker{ii}='B+';
                case 'o',           marker{ii}='Bo';
                case '*',           marker{ii}='Basterisk';
                case '.',           marker{ii}='*';
                case 'x',           marker{ii}='B|';        %%%%%% different
                case 'square',      marker{ii}='Bsquare';
                case 'diamond',     marker{ii}='diamond';
                case '^',           marker{ii}='Btriangle';
                case 'v',           marker{ii}='square*';   %%%%%% different
                case '>',           marker{ii}='diamond*';  %%%%%% different
                case '<',           marker{ii}='triangle*'; %%%%%% different
                case 'pentagram',   marker{ii}='Bpentagon';
                case 'hexagram',    marker{ii}='pentagon*'; %%%%%% different
                otherwise
                    marker{ii}='none';
            end
        end
    end


    function [fig,globals] = update_figure_geometry(fig,globals)
        % This function generstes a gui, where some parameters can be set, the
        % figure is redrawn afterwards.
        
        exit_flag=0; % init
        % AXES:
        set(fig.axh,'Units','pixels'); ap=get(fig.axh,'Position');
        if ~isempty(ahSecond)
            set(fig.axh2,'Units','pixels'); ap2=get(fig.axh2,'Position');
        end
        
        % FIGURE:
        set(fig.figh,'Units','pixels');
        fp=get(fig.figh,'Position');
        topSpace=fp(4)-ap(4)-ap(2);
        globals.FigAxWidthFac=fp(3)/ap(3);
        
        % Create GUI figure
%         guisize=[260 fp(4)]; % pixel (300 x 200)
        guisize=[310 420];
        guiPos=[fp(1)+fp(3),fp(2),guisize];
        
        if isempty(globals.heightFactor), % set height factor automatically if not specified
            globals.heightFactor=ap(4)/ap(3);
        else
            update; % update figure
        end
        
        fGui = figure('units','pixels','Position',guiPos,...
            'menubar','none','name','Figure Export Properties',...
            'numbertitle','off','resize','off','WindowStyle','modal',...
            'Color',get(0,'defaultuicontrolbackgroundcolor'),...
            'CloseRequestFcn',@UImenuExit_Callback);
        % panels
        uipanelGeometry = uipanel('parent',fGui,'units','pixels',...
            'position',[10 150 290 260],'title','Geometry Settings',...
            'FontSize',9,'FontWeight','bold');
        uipanelAuto = uipanel('parent',uipanelGeometry,'units','pixels',...
            'position',[10 10 130 130],'title','Auto (set by textsize)',...
            'FontSize',8,'FontWeight','normal');
        uipanelManual = uipanel('parent',uipanelGeometry,'units','pixels',...
            'position',[150 10 130 130],'title','Manual',...
            'FontSize',8,'FontWeight','normal');
        uipanelMisc = uipanel('parent',fGui,'units','pixels',...
            'position',[10 10 290 130],'title','Misc',...
            'FontSize',9,'FontWeight','bold');
        
        % fields
        WidthTextPre = uicontrol(uipanelGeometry,'style','text','units','pixels',...
            'position',[10 210 75 17],'string','Width / cm:',...
            'HorizontalAlignment','left');
        editWidth = uicontrol(uipanelGeometry,'style','edit','units','pixels',...
            'position',[90 210 50 20],'string',globals.plotwidth,'BackgroundColor','white',...
            'Enable','on',...
            'Callback',@editWidth_callback);
        if strcmp(globals.widthType,'axes'), % test if axes is activated
            axtrue=1; figtrue=0; else axtrue=0; figtrue=1;
        end
        radiobuttonFig = uicontrol(uipanelGeometry,'style','radiobutton',...
            'units','pixels','position',[150 210 50 20],...
            'string','Figure','Value',figtrue,'Enable','on',...
            'Callback',@radiobuttonFig_callback);
        radiobuttonAxes = uicontrol(uipanelGeometry,'style','radiobutton',...
            'units','pixels','position',[225 210 50 20],...
            'string','Axes','Value',axtrue,'Enable','on',...
            'Callback',@radiobuttonAxes_callback);
        textFactor = uicontrol(uipanelGeometry,'style','text',...
            'units','pixels','position',[10 180 75 17],...
            'string','Height /  Width:','HorizontalAlignment','left');
        editFactor = uicontrol(uipanelGeometry,'style','edit','units','pixels',...
            'position',[90 180 50 20],'string',globals.heightFactor,'BackgroundColor','white',...
            'Enable','on',...
            'Callback',@editFactor_callback);
        radiobuttonAspectRatio = uicontrol(uipanelGeometry,'style','radiobutton',...
            'units','pixels','position',[150 180 80 20],...
            'string','Axis equal','Value',0,...
            'Callback',@radiobuttonAspectRatio_callback);
        radiobuttonSquare = uicontrol(uipanelGeometry,'style','radiobutton',...
            'units','pixels','position',[225 180 55 20],...
            'string','Square','Value',0,...
            'Callback',@radiobuttonSquare_callback);
        textSpaceSurrounding = uicontrol(uipanelGeometry,'style','text',...
            'units','pixels','position',[10 145 180 17],'FontWeight','bold',...
            'string','Space surrounding the axes:','HorizontalAlignment','left');
        
        sliderFont = uicontrol(uipanelAuto,'style','slider',...
            'units','pixels','position',[10 85 105 20],...
            'string','Figure','BackgroundColor','white',...
            'Min',globals.labelsepMin,'Max',globals.labelsepMax,...
            'Value',globals.labelsep,'SliderStep',[1/(globals.labelsepMax-globals.labelsepMin),0],...
            'Callback',@sliderFont_callback);
        radiobuttonCenterAxes = uicontrol(uipanelAuto,'style','radiobutton',...
            'units','pixels','position',[10 60 80 20],...
            'string','Center axes','Value',globals.centerAxes,'Enable','on',...
            'Callback',@radiobuttonCenterAxes_callback);
        texthorizontally = uicontrol(uipanelAuto,'style','text',...
            'units','pixels','position',[28 45 75 17],...
            'string','horizontally','HorizontalAlignment','left');
        
        LeftPre = uicontrol(uipanelManual,'style','text','units','pixels',...
            'position',[10 85 37 17],'string','Left:','HorizontalAlignment','left');
        BottomPre = uicontrol(uipanelManual,'style','text','units','pixels',...
            'position',[10 60 37 17],'string','Bottom:','HorizontalAlignment','left');
        RightPre = uicontrol(uipanelManual,'style','text','units','pixels',...
            'position',[10 35 37 17],'string','Right:','HorizontalAlignment','left');
        TopPre = uicontrol(uipanelManual,'style','text','units','pixels',...
            'position',[10 10 37 17],'string','Top:','HorizontalAlignment','left');
        LeftPost = uicontrol(uipanelManual,'style','text','units','pixels',...
            'position',[105 85 20 17],'string','cm','HorizontalAlignment','left');
        BottomPost = uicontrol(uipanelManual,'style','text','units','pixels',...
            'position',[105 60 20 17],'string','cm','HorizontalAlignment','left');
        RightPost = uicontrol(uipanelManual,'style','text','units','pixels',...
            'position',[105 35 20 17],'string','cm','HorizontalAlignment','left');
        TopPost = uicontrol(uipanelManual,'style','text','units','pixels',...
            'position',[105 10 20 17],'string','cm','HorizontalAlignment','left');
        globals.ManualTightStrings = {'---','---','---','---'}; % init
        globals.tightInset_cmBackup= globals.tightInset_cm; % copy
        if ~ isempty (globals.tightInset_cm)
            for ii=1:4,
                globals.ManualTightStrings{ii} = num2str(globals.tightInset_cm (ii));
            end
        else
            globals.tightInset_cmBackup = [1 1 1 1];
        end
        if axtrue, EnableTight = 'on'; else EnableTight = 'off'; end
        editLeft = uicontrol(uipanelManual,'style','edit','units','pixels',...
            'position',[50 85 50 20],'string',globals.ManualTightStrings{1},'BackgroundColor','white',...
            'Enable','on','Enable',EnableTight,'Callback',@editLeft_callback);
        editBottom = uicontrol(uipanelManual,'style','edit','units','pixels',...
            'position',[50 60 50 20],'string',globals.ManualTightStrings{2},'BackgroundColor','white',...
            'Enable','on','Enable',EnableTight,'Callback',@editBottom_callback);
        editRight = uicontrol(uipanelManual,'style','edit','units','pixels',...
            'position',[50 35 50 20],'string',globals.ManualTightStrings{3},'BackgroundColor','white',...
            'Enable','on','Enable',EnableTight,'Callback',@editRight_callback);
        editTop = uicontrol(uipanelManual,'style','edit','units','pixels',...
            'position',[50 10 50 20],'string',globals.ManualTightStrings{4},'BackgroundColor','white',...
            'Enable','on','Enable',EnableTight,'Callback',@editTop_callback);
        
        radiobuttonBoundingBox = uicontrol(uipanelMisc,'style','radiobutton',...
            'units','pixels','position',[10 90 120 20],...
            'string','Draw figure outline','Value',globals.BoundingBox,...
            'Callback',@radiobuttonBoundingBox_callback);
        radiobuttonStandAloneFile = uicontrol(uipanelMisc,'style','radiobutton',...
            'units','pixels','position',[170 90 100 20],...
            'string','Stand-alone file','Value',globals.StandAloneFile,...
            'Callback',@radiobuttonStandAloneFile_callback);
        radiobuttonReducePoints = uicontrol(uipanelMisc,'style','radiobutton',...
            'units','pixels','position',[10 68 170 20],...
            'string','Reduce data points, Factor:','Value',1,...
            'Callback',@radiobuttonReducePoints_callback);
        editReducePointsFactor = uicontrol(uipanelMisc,'style','edit','units','pixels',...
            'position',[170 68 50 20],'string',globals.PointReduceFactor,'BackgroundColor','white',...
            'Enable','on',...
            'Callback',@editReducePointsFactor_callback);
        LegendText = uicontrol(uipanelMisc,'style','text','units','pixels',...
            'position',[10 40 75 17],'string','Legend style:',...
            'HorizontalAlignment','left');
        
        
        globals.LegendStyle = {'Rectangular Box', ...
            'Rounded corners', 'Shadowbox', 'No Box'};
        globals.LegendBoxTeXcode = get_LegendTexCode_Pre (globals.idx_LegendStyle);
        
        popupLegendStyle = uicontrol(uipanelMisc,'Style','popupmenu',...
            'units','pixels','position',[80 40 120 20],...
            'String', globals.LegendStyle,'BackgroundColor','white',...
            'Value', globals.idx_LegendStyle, ...
            'Callback',@popupLegendStyle_callback);
        SaveText = uicontrol(uipanelMisc,'style','text','units','pixels',...
            'position',[10 10 75 17],'string','Export to file:',...
            'HorizontalAlignment','left');
        popupSaveFile = uicontrol(uipanelMisc,'Style','popupmenu',...
            'units','pixels','position',[80 10 120 20],...
            'String', {globals.filename,'Select...'},'BackgroundColor','white',...
            'Callback',@popupSaveFile_callback);
        pushbuttonDone = uicontrol(uipanelMisc,'style','pushbutton','units','pixels',...
            'position',[210 10 70 30],'string','Start',...
            'Callback',@DoneButton_callback);
        
        if globals.ShowGUI
            uiwait; % wait until figure is deleted
        else
            DoneButton_callback;
        end
        
        function editWidth_callback(~,~)
            globals.plotwidth=str2num(get(editWidth,'string'));
        end
        
        function radiobuttonAxes_callback(~,~)
            isOn=get(radiobuttonAxes,'Value');
            if isOn,
                set(radiobuttonFig,'Value',0)
                globals.widthType='axes';
                set(editLeft,'Enable','on');  set(editBottom,'Enable','on');
                set(editRight,'Enable','on'); set(editTop,'Enable','on');
            else
                set(radiobuttonFig,'Value',1)
                globals.widthType='figure';
                set(editLeft,'Enable','off');  set(editBottom,'Enable','off');
                set(editRight,'Enable','off'); set(editTop,'Enable','off');
            end
        end
        
        function radiobuttonFig_callback(~,~)
            isOn=get(radiobuttonFig,'Value');
            if isOn,
                set(radiobuttonAxes,'Value',0)
                globals.widthType='figure';
                set(editLeft,'Enable','off');  set(editBottom,'Enable','off');
                set(editRight,'Enable','off'); set(editTop,'Enable','off');
            else
                set(radiobuttonAxes,'Value',1)
                globals.widthType='axes';
                set(editLeft,'Enable','on');  set(editBottom,'Enable','on');
                set(editRight,'Enable','on'); set(editTop,'Enable','on');
            end
        end
        
        function radiobuttonCenterAxes_callback(~,~)
            isOn=get(radiobuttonCenterAxes,'Value');
            if isOn,
                globals.centerAxes=1;
            else
                globals.centerAxes=0;
            end
            globals.tightInset_cm = []; % backup still exist
            set(editLeft,'String',globals.ManualTightStrings{1});
            set(editBottom,'String',globals.ManualTightStrings{2});
            set(editRight,'String',globals.ManualTightStrings{3});
            set(editTop,'String',globals.ManualTightStrings{4});
        end
        
        function radiobuttonReducePoints_callback(~,~)
            reduceOn=get(radiobuttonReducePoints,'Value');
            if reduceOn,
                set(radiobuttonReducePoints,'Value',1)
                set(editReducePointsFactor,'Enable','on')
                globals.PointReduceFactor=str2num(get(editReducePointsFactor,'string'));
            else
                set(radiobuttonReducePoints,'Value',0)
                set(editReducePointsFactor,'Enable','off')
                globals.PointReduceFactor=0; % reduce off
            end
        end
        
        function editReducePointsFactor_callback(hObject,eventdata)
            globals.PointReduceFactor=str2num(get(editReducePointsFactor,'string'));
        end
        
        function editFactor_callback(hObject,eventdata)
            globals.heightFactor=str2num(get(editFactor,'string'));
            set(radiobuttonAspectRatio,'Value',0)
            set(radiobuttonSquare,'Value',0)
            update;
        end
        
        function radiobuttonAspectRatio_callback(hObject,eventdata)
            globals.heightFactor=(fig.YLim(2)-fig.YLim(1))/(fig.XLim(2)-fig.XLim(1));
            set(editFactor,'string',globals.heightFactor);
            set(radiobuttonSquare,'Value',0)
            update;
        end
        
        function radiobuttonSquare_callback(hObject,eventdata)
            globals.heightFactor=1;
            set(editFactor,'string',globals.heightFactor);
            set(radiobuttonAspectRatio,'Value',0)
            update;
        end
        
        function editLeft_callback(hObject,eventdata)
            globals.tightInset_cm(1) = str2num(get(editLeft,'string'));
            globals.tightInset_cmBackup(1) = globals.tightInset_cm(1);
            globals.tightInset_cm(2:4) = globals.tightInset_cmBackup(2:4);
            globals.centerAxes = 0; set(radiobuttonCenterAxes,'Value',0); % deactivate
            set(editLeft,'string',globals.tightInset_cm(1)); % update
            set(editBottom,'string',globals.tightInset_cm(2)); % update
            set(editRight,'string',globals.tightInset_cm(3)); % update
            set(editTop,'string',globals.tightInset_cm(4)); % update
            update;
        end
        
        function editBottom_callback(hObject,eventdata)
            globals.tightInset_cm(2) = str2num(get(editBottom,'string'));
            globals.tightInset_cmBackup(2) = globals.tightInset_cm(2);
            globals.tightInset_cm(1) = globals.tightInset_cmBackup(1);
            globals.tightInset_cm(3:4) = globals.tightInset_cmBackup(3:4);
            globals.centerAxes = 0; set(radiobuttonCenterAxes,'Value',0); % deactivate
            set(editLeft,'string',globals.tightInset_cm(1)); % update
            set(editBottom,'string',globals.tightInset_cm(2)); % update
            set(editRight,'string',globals.tightInset_cm(3)); % update
            set(editTop,'string',globals.tightInset_cm(4)); % update
            update;
        end
        
        function editRight_callback(~,eventdata)
            globals.tightInset_cm(3) = str2num(get(editRight,'string'));
            globals.tightInset_cmBackup(3) = globals.tightInset_cm(3);
            globals.tightInset_cm(1:2) = globals.tightInset_cmBackup(1:2);
            globals.tightInset_cm(4) = globals.tightInset_cmBackup(4);
            globals.centerAxes = 0; set(radiobuttonCenterAxes,'Value',0); % deactivate
            set(editLeft,'string',globals.tightInset_cm(1)); % update
            set(editBottom,'string',globals.tightInset_cm(2)); % update
            set(editRight,'string',globals.tightInset_cm(3)); % update
            set(editTop,'string',globals.tightInset_cm(4)); % update
            update;
        end
        
        function editTop_callback(~,~)
            globals.tightInset_cm(4) = str2num(get(editTop,'string'));
            globals.tightInset_cmBackup(4) = globals.tightInset_cm(4);
            globals.tightInset_cm(1:3) = globals.tightInset_cmBackup(1:3);
            globals.centerAxes = 0; set(radiobuttonCenterAxes,'Value',0); % deactivate
            set(editLeft,'string',globals.tightInset_cm(1)); % update
            set(editBottom,'string',globals.tightInset_cm(2)); % update
            set(editRight,'string',globals.tightInset_cm(3)); % update
            set(editTop,'string',globals.tightInset_cm(4)); % update
            update;
        end
        
        function sliderFont_callback(~,~)
            globals.labelsep=get(sliderFont,'Value'); % value
            set(fig.axh,'FontSize',globals.labelsep);
            set(get(fig.axh,'Title'),'FontSize',globals.labelsep);
            set(get(fig.axh,'Xlabel'),'FontSize',globals.labelsep);
            set(get(fig.axh,'Ylabel'),'FontSize',globals.labelsep);
            if ~isempty(ahSecond)
                set(fig.axh2,'FontSize',globals.labelsep);
                set(get(fig.axh2,'Ylabel'),'FontSize',globals.labelsep);
            end
            globals.tightInset_cm = []; % backup still exist
            set(editLeft,'String',globals.ManualTightStrings{1});
            set(editBottom,'String',globals.ManualTightStrings{2});
            set(editRight,'String',globals.ManualTightStrings{3});
            set(editTop,'String',globals.ManualTightStrings{4});
            update;
        end
        
        function radiobuttonBoundingBox_callback(~,~)
            if get(radiobuttonBoundingBox,'Value') == 1, %
                globals.BoundingBox=1;
            else
                globals.BoundingBox=0;
            end
        end
        
        function radiobuttonStandAloneFile_callback(~,~)
            if get(radiobuttonStandAloneFile,'Value') == 1, %
                globals.StandAloneFile=1;
            else
                globals.StandAloneFile=0;
            end
        end
        
        
        
        function out = get_LegendTexCode_Pre (val)
            switch val
                case 1, % rectangular box
                    out = {'\psframebox[framesep=0pt,linewidth=\AxesLineWidth]{\psframebox*{\begin{tabular}{l}'}; % init text for framebox
                case 2, % rounded corners
                    out = {'\psframebox[fillstyle=solid,cornersize=absolute,linearc=5pt,framesep=0pt,linewidth=\AxesLineWidth]{\psframebox[linestyle=none]{\begin{tabular}{l}'}; % init text for rounded corners
                case 3, % shadowbox
                    out = {'\psshadowbox[framesep=0pt,linewidth=\AxesLineWidth]{\psframebox*{\begin{tabular}{l}'}; % init text for shadowbox
                case 4, % no box
                    out = {'{\psframebox*[framesep=0pt]{\begin{tabular}{l}'}; % init text for framebox
            end
        end
        
        function popupLegendStyle_callback(~,~)
            globals.idx_LegendStyle = get(popupLegendStyle,'Value');
            globals.LegendBoxTeXcode = get_LegendTexCode_Pre (globals.idx_LegendStyle);
        end
        
        function popupSaveFile_callback(~,~)
            if get(popupSaveFile,'Value')==2, % choose new filename
                [filename,pathname] = uiputfile(fullfile(globals.pathname,globals.filename));
                if ischar(filename) , % correct file
                    globals.filename=filename;
                    globals.pathname=pathname;
                end
                set(popupSaveFile,'String',{globals.filename,'Select...'})
                set(popupSaveFile,'Value',1)
            end
        end
        
        function DoneButton_callback(~,~)
            if get(radiobuttonAxes,'Value')==1, % width for axes
                globals.FigAxWidthFac=1;
            end
            delete(fGui);  % delete figure and exit
        end
        
        function UImenuExit_Callback(~,~)
            close(fig.figh)
            delete(fGui)
            exit_flag=1;
        end
        
        function update % update figure
            % update figure and axes
            if not(isempty(fig.LegendLoc)),
                % if a legend is present, it must first be disabled when
                % changing the figure size, it is replotted after change in
                % figure size
                legend(fig.axh,'off')
            end
            ap(4)=globals.heightFactor*ap(3); % update height of axes
            fp(4)=topSpace+ap(4)+ap(2); % update height of figure
            set(fig.figh,'Units','pixels','Position',round(fp));
            set(fig.axh,'Units','pixels','Position',round(ap));
            set(fig.axh,'XTickLabel',fig.XTickLabel);
            xlabel(fig.axh,[fig.XLabel,fig.XLabelPost],'FontSize',globals.labelsep,'interpreter','none');
            ylabel(fig.axh,[fig.YLabel,fig.YLabelPost],'FontSize',globals.labelsep,'interpreter','none');
            title(fig.axh,fig.Title,'FontSize',globals.labelsep,'interpreter','none');
            if not(isempty(fig.LegendLoc)) && isempty(regexp(fig.LegendLoc,'Outside', 'once')), % a legend is replotted if the legend is not outside
                legend(fig.axh,fig.LegendStr,'Location',fig.LegendLoc,'interpreter','none')
            end
            if ~isempty(ahSecond)
                ap2(4)=globals.heightFactor*ap2(3); % update height of axes
                set(fig.axh2,'Position',ap2);
                ylabel(fig.axh2,[fig.YLabel2,fig.YLabelPost2],'FontSize',globals.labelsep,'interpreter','none');
            end
            % movegui(fig.figh,'center'); % shift figure to screen center
        end
        
        if not(exit_flag)
            update
            set(fig.axh,'Units','normalized');
        end
    end


% Copied from Ercan Solak's file fig2tex.m
    function resstr = strfill(genstr,fpar)
        %STRFILL Replace the numbered tokens with parameters
        %   STRFILL(GENSTR,FPAR) replaces the numbered token
        %   #i# in the string GENSTR with the ith element
        %   of the cellarray FPAR. This script is used by
        %   FIG2TEX.
        
        resstr=genstr;
        for ii_str=1:length(fpar)
            if isnumeric(fpar{ii_str})
                reptoken = num2str(fpar{ii_str},globals.FloatingPointFormat);
            else
                reptoken = fpar{ii_str};
            end;
            resstr = strrep(resstr,['#',num2str(ii_str),'#'],reptoken);
        end;
        
    end


% FUNCTION REDUCE_LINE
%
%  COPY OF MATFIG2PGF

%--------------------------------------------------------------------------
    function [ reduced_xydata ] = reduce_line( xydata, maxerror )
        
        
        N = size(xydata,2);
        
        if (maxerror <= 0) || (N < 3)
            reduced_xydata = xydata;
            return
        end
        
        xydata = minmaxdecimation(xydata, 4*maxerror);
        
        N = size(xydata,2);
        
        plotPoints = 1;
        lastPlotted = 1;
        
        xdata = xydata(1,:);
        ydata = xydata(2,:);
        
        for i = 3:N
            % Calculate distance
            % see http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
            % x1 = xydata(:,i)  x2 = xydata(:,lastPlotted)  x0 = all points
            % between x1 and x2
            % See also: http://astronomy.swin.edu.au/~pbourke/geometry/pointline/
            p1 = xydata(:,i);
            p2 = xydata(:,lastPlotted);
            dp = sqrt(sum((p1-p2).^2));
            frac1 = ( xdata(lastPlotted+1:i-1)-xdata(i) ) .* ( xdata(lastPlotted)-xdata(i) );
            frac1 = frac1 + ( ydata(lastPlotted+1:i-1)-ydata(i) ) .* ( ydata(lastPlotted)-ydata(i) );
            u = frac1 ./ sum((p2-p1).^2);
            
            % Convert u to the distance from the point to p1 or p2
            % For points where the closest point on line p1-p2 is outside of p1 and
            % p2 u is now the distance to p1 or p2. When the closest point was
            % between p1 and p2 u will be zero.
            u((u > 0) & (u <= 1)) = 0;
            u(u>1) = dp*(u(u>1)-1);
            u(u<0) = dp*(-u(u<0));
            
            % Calculate shortest distance from point to line (p1-p2)
            a = xdata(lastPlotted)-xdata(i);
            b = ydata(lastPlotted)-ydata(i);
            c = xdata(i)-xdata(lastPlotted+1:i-1);
            d = ydata(i)-ydata(lastPlotted+1:i-1);
            frac1 = abs(a.*d-b.*c);
            frac2 = sqrt(sum((xydata(:,lastPlotted)-xydata(:,i)).^2));
            d = frac1./frac2;
            
            d = sqrt(d.^2 + u.^2);
            
            if max(d) > maxerror
                lastPlotted = i-1;
                plotPoints = [plotPoints lastPlotted];
            end
            
        end
        
        plotPoints = [plotPoints N];
        if N > 5
            reduced_xydata = xydata(:,plotPoints);
        else
            reduced_xydata = xydata;
        end
        if N ~= size(reduced_xydata,2)
            fprintf('reduced data points: org: %d -> red: %d\n',...
                N, size(reduced_xydata,2));
        end
        
        %-- end of function reduce_line -------------------------------------------
        
    end

% FUNCTION MINMAXDECIMATION
%
%  COPY OF MATFIG2PGF
%
    function [ decimatedXydata ] = minmaxdecimation( xydata, columnWidth )
        
        xdata = xydata(1,:);
        ydata = xydata(2,:);
        minx = min(xdata);
        maxx = max(xdata);
        N = (maxx-minx)/columnWidth;  % number of columns
        
        % dx = xdata(i+1)-xdata(i) is the same for all values. Otherwise this
        % function is not going to work.
        dx = diff(xdata);
        maxdx = max(dx);
        mindx = min(dx);
        dx = median(dx);
        
        % If dx is not the same for all values OR
        % if the number of columns is less than 0.5*data length
        % then we can not use min-max decimation
        if ( (maxdx-mindx)/dx > 1e-6 ) || ( N > 0.5*length(xdata) )
            decimatedXydata = xydata;
            return;
        end
        
        decimatedXydata = [];
        lastIndex = 0;
        for i = 1:ceil(N)
            thisIndex = min( floor(i*columnWidth/dx), length(xdata) );
            [miny, minIndex] = min(ydata(lastIndex+1:thisIndex));
            [maxy, maxIndex] = max(ydata(lastIndex+1:thisIndex));
            minIndex = minIndex+lastIndex;
            maxIndex = maxIndex+lastIndex;
            if minIndex < maxIndex
                decimatedXydata = [decimatedXydata [xdata(minIndex);miny] [xdata(maxIndex);maxy]];
            else
                decimatedXydata = [decimatedXydata [xdata(maxIndex);maxy] [xdata(minIndex);miny]];
            end
            lastIndex = thisIndex;
        end
        if size(decimatedXydata,2) > 10
            fprintf('  min-max decimation: original %d  decimated %d\n', size(xydata,2), size(decimatedXydata,2));
        else
            decimatedXydata = xydata;
        end
        %- end of function minmaxdecimation ---------------------------------------
    end




    function [Pcut, Psnip] = cut_data (xd, yd, XLim, YLim)
        % cut data that are out of bounds and append them as extra lines:
        
        % Pcut is a cell array with each element being a nx2 matrix
        % containing a point in each row that belongs to a valid line
        % segment
        %
        % Psnip is a cell array with each element being a nx2 matrix
        % containing a point in each row that belongs to a valid snipped
        % that should be appended at the beginning or end
        
        % init
        P_tmp = [];
        Pcut  = {};
        Psnip = {}; % here come all snippet points
        
        % scan through all points
        for ii = 2 : length(xd), % start with second point
            
            P_curr = [xd(ii),   yd(ii)];
            P_last = [xd(ii-1), yd(ii-1)];
            
            if ~isfinite(P_curr(1)) || ~isfinite(P_curr(2)), % at least one value is nan
                
                if ~isempty(P_tmp),
                    Pcut{end+1} = P_tmp; % add last segments to cell, ...
                    P_tmp = []; % and reset
                end                
                
            else
                                
                if ~P_is_in(P_curr), % not in
                    
                    if P_is_in(P_last) && ... % last point was in, new snippet
                            (isfinite(P_last(1)) && isfinite(P_last(2))), % valid
                        
                        Psnip{end+1} = [P_last; interpolate(P_curr, P_last)];
                    elseif ~((P_is_above(P_curr) && P_is_above(P_last)) || ...  % both P_curr and P_last are out, but are above and below
                            (P_is_below(P_last) && P_is_below(P_curr))) 
                        Psnip{end+1} = interpolate(P_curr, P_last);      
                        Psnip{end}
                    end
                                                                                
                    if ~isempty(P_tmp),
                        Pcut{end+1} = P_tmp; % add last segments to cell, ...
                        P_tmp = []; % and reset
                    end
                    
                elseif P_is_in(P_curr),
                    
                    P_tmp(end+1,:) = P_curr; % add this to main points
                    
                    if ~P_is_in(P_last) && ... % last point was not in, new snippet
                            (isfinite(P_last(1)) && isfinite(P_last(2))), % valid
                        
                        Psnip{end+1} = [interpolate(P_curr, P_last); P_curr];
                        
                    end
                    
                    % treat first and last segments special
                    
                    if ii == 2, % first element
                    
                        if P_is_in(P_last) && ... % last point was in append this point:
                            (isfinite(P_last(1)) && isfinite(P_last(2))), % valid
                            
                            P_tmp = [P_last; P_tmp];
                            
                        end
                        
                    end
                    
                    if ii == length(xd), % last element, write data
                        
                        Pcut{end+1} = P_tmp;
                        
                    end
                    
                    
                end
            end
            
        end
        
        
        % Subfunctions:
        
        function [Pint, varargout] = interpolate(P1, P2)
        
            cut_pos = ''; % init, left, right, top, bottom
            
            if P_is_in(P1) && P_is_in(P2)
                Pint = []; % both in, nothing to interpolate
            elseif ~P_is_in(P1) && ~P_is_in(P2) % both out, return both points                
                % build straight line (y = m * x + b)
                m = (P2(2) - P1(2)) / (P2(1) - P1(1));
                b = P2(2) - m * P2(1);
                % (x,YLim(1))
                dx1 = (YLim(1) - b) / m; 
                if dx1 < XLim(1), dx1 = XLim(1);
                elseif dx1 > XLim(2), dx1 = XLim(2); end
                
                % (x,YLim(2))
                dx2 = (YLim(2) - b) / m;
                if dx2 < XLim(1), dx2 = XLim(1);
                elseif dx2 > XLim(2), dx2 = XLim(2); end
                Pint = [dx1,YLim(1);dx2,YLim(2)];
                cut_pos = 'top-bottom';
            else
 
                % build straight line (y = m * x + b)
                m = (P2(2) - P1(2)) / (P2(1) - P1(1));
                b = P2(2) - m * P2(1);
                
                % test where straght line intersects with borders
                if ((P1(1) >= XLim(1)) && (P2(1) < XLim(1))) || ...
                      ((P2(1) >= XLim(1)) && (P1(1) < XLim(1))), % left
                  % 'left'
                  Pint = [XLim(1), m*XLim(1)+b];
                  cut_pos = 'left';
                  
                elseif ((P1(1) > XLim(2)) && (P2(1) <= XLim(2))) || ...
                      ((P2(1) > XLim(2)) && (P1(1) <= XLim(2))), % left
                  % 'right'
                  Pint = [XLim(2), m*XLim(2)+b];
                  cut_pos = 'right';
                  
                elseif ((P1(2) >= YLim(1)) && (P2(2) < YLim(1))) || ...
                        ((P2(2) >= YLim(1)) && (P1(2) < YLim(1))), % left
                    % 'bottom'
                    % test for vertical line:
                    if ~isfinite(m), % vertical
                        Pint = [P2(1), YLim(1)];
                    else
                        Pint = [(YLim(1)-b)/m, YLim(1)];
                    end
                    cut_pos = 'bottom';

                elseif ((P1(2) > YLim(2)) && (P2(2) <= YLim(2))) || ...
                      ((P2(2) > YLim(2)) && (P1(2) <= YLim(2))), % left
                  % 'top'
                  if ~isfinite(m), % vertical
                      Pint = [P2(1), YLim(2)];
                  else
                      Pint = [(YLim(2)-b)/m, YLim(2)];
                  end
                  cut_pos = 'top';
                  
                else
                    
                    warning('Something went wrong')
                    
                end
               
            end
            
            nout = max(nargout,1)-1;
            if nout == 1,
                varargout(1) = {cut_pos};
            end

            
        end
        
        function answer = P_is_in(P)
            answer = false;
            x = P(1);
            y = P(2);
            if (x >= XLim(1)) && (x <= XLim(2)) && ...
                    (y >= YLim(1)) && (y <= YLim(2))
                answer = true;
            end
        end
        
        function answer = P_is_above(P)
            answer = false;
            x = P(1);
            y = P(2);
            if (x <= XLim(1) || y >= YLim(1))
                answer = true;
            end
        end
        
        function answer = P_is_below(P)
            answer = false;
            x = P(1);
            y = P(2);
            if (x >= XLim(2) || y <= YLim(2))
                answer = true;
            end
        end
        
        
        
    end




% FUNCTION TEXT_CONTAINS_DOLLARSYMBOL
%
%    answer = text_contains_dollarsymbol(textin)
%
%       textin  - normal text to be tested for at least one '$'
%
%       answer - boolean answer
%--------------------------------------------------------------------------
    function answer = text_contains_dollarsymbol(textin)
        answer=false;
        if ~iscell(textin), textin={textin}; end % convert to cell
        for ii=1:length(textin),
            if iscell(textin{ii}) % if still cell, multiple lines
                for jj=1:length(textin{ii}), % another loop for two line entries
                    if ~isempty(strfind(textin{ii}{jj},'$')), % at least one $ has been found
                        answer=true;
                        break
                    end
                end
            else % just one line
                if ~isempty(strfind(textin{ii},'$')), % at least one $ has been found
                    answer=true;
                    break
                end
            end
            
        end
    end

% FUNCTION CONVERT_TEXT_TO_LATEX
%
%    textout = convert_text_to_line(textin)
%
%       textin  - normal text to be converted to proper latex code
%
%       textout - converted text
%
%  COPY OF MATFIG2PGF
%
%--------------------------------------------------------------------------
    function textout = convert_text_to_latex(textin)
        textcell=textin;
        % test if textin is a cell array
        if ~iscell(textin), textcell={textin}; end % convert
        textoutCell=textcell;
        for jj=1:length(textcell),
            textin=textcell{jj};
            if iscell(textin), % still a cell, means multiple lines
                for kk=1:length(textin)
                    textout{kk}=doTheConversion(textin{kk});
                end
                textoutCell{jj}=textout;
            else
                textout=doTheConversion(textin);
                textoutCell{jj}=textout;
            end
        end
        
        if length(textcell)==1,
            textout=textoutCell{1};
        else
            textout=textoutCell;
        end
        
        %-- end of function convert_text_to_latex ---------------------------------
        
        
        function textout=doTheConversion(textin)
            % textin
            
            % Split strings in groups separated by space and/or }[a-z]
            splitStrings = {};
            i = 1;
            thisStartIndex = 1;
            while i <= length(textin),
                if ~isempty(regexp(textin(i), '\s', 'once'))
                    splitStrings{length(splitStrings)+1} = textin(thisStartIndex:i);
                    thisStartIndex = i+1;
                elseif (i < length(textin)) && ~isempty(regexpi(textin(i:i+1), '}[a-z]', 'once'))
                    splitStrings{length(splitStrings)+1} = textin(thisStartIndex:i);
                    thisStartIndex = i+1;
                elseif (i < length(textin)) && ~isempty(regexpi(textin(i:i+1), '[^_\^]{'))
                    splitStrings{length(splitStrings)+1} = textin(thisStartIndex:i);
                    thisStartIndex = i+1;
                elseif i == length(textin)
                    % Last character of string
                    splitStrings{length(splitStrings)+1} = textin(thisStartIndex:i);
                end
                i = i+1;
            end
            
            % If two consecutive strings need to set in mathmode and have no whitespace
            % in between. They must be joined to one math mode string.
            newSplitStrings = {};
            for i = 1:length(splitStrings)
                if i > 1
                    prev = newSplitStrings{length(newSplitStrings)};
                    next = splitStrings{i};
                    if inMathMode(prev) && inMathMode(next)
                        if isempty(regexp(prev(end), '\s', 'once')) && isempty(regexp(next(1), '\s', 'once'))
                            newSplitStrings{length(newSplitStrings)} = [prev next];
                        else
                            newSplitStrings{length(newSplitStrings)+1} = next;
                        end
                    else
                        newSplitStrings{length(newSplitStrings)+1} = next;
                    end
                else
                    newSplitStrings{length(newSplitStrings)+1} = splitStrings{i};
                end
            end
            splitStrings = newSplitStrings;
            
            textout = '';
            for i = 1:length(splitStrings)
                if iscell(splitStrings{i}), % not sure if this is the right way
                    splitStrings{i}=splitStrings{i}{1};
                end
                if ~inMathMode(splitStrings{i})
                    textout = [textout splitStrings{i}];
                else
                    thisString = splitStrings{i};
                    
                    % Remove whitespace at end of string
                    lastIndex = length(thisString)+1-regexp(fliplr(thisString), '[^\s]', 'once');
                    if lastIndex < length(thisString)
                        trailingWhitespace = true;
                        thisString = thisString(1:lastIndex);
                    else
                        trailingWhitespace = false;
                    end
                    
                    % If the are acculades at the beginning and end they can be removed
                    if strcmp(thisString(1), '{') && strcmp(thisString(lastIndex), '}')
                        thisString = thisString(2:lastIndex-1);
                    end
                    
                    textout = [textout '$' thisString '$'];
                    if trailingWhitespace
                        textout = [textout ' '];
                    end
                end
            end
            
            % Replace % signs in the text because they are comments in latex
            textout = regexprep(textout, '%', '\\%');
        end
        
        
    end

% FUNCTION INMATHMODE
%
% Determines whether a string needs to be typeset in math mode in LaTeX
%
% [ mathmode ] = inMathMode( str )
%
%  str - string that needs to be checked
%
%  mathmode - True when it needs to be typeset in math mode, return false
%             when it should be typeset in normal text mode.
%
%  COPY OF MATFIG2PGF
%
    function [ mathmode ] = inMathMode( str )
        mathmode = ~isempty(regexp(str, '[|\\_\^]', 'once'));
        %-- end of function inmathmode --------------------------------------------
    end





    function labels_out = convert_TickLabels2Cell (ticks, labels)
        labels_out = {};
        if length (ticks) == size (labels, 1), % same size
            
            if ~ iscell (labels)
                for ii = 1 : size (labels, 1)
                    labels_out{ii} = labels(ii, :);
                end
            else
                labels_out = labels;
            end
            
        else % not the same size
            
            if ~isempty(ticks) && isempty(labels), % ticks but no labels,
                
                for ii = 1 : length (ticks)
                    labels_out{ii} = '';% empty label
                end
                
            else % not same size but 
   
                reps = ceil (length (ticks) / size (labels, 1));
                
                if ~ iscell (labels)
                    
                    labels = repmat (labels, reps, 1); % repeat labels
                    for ii = 1 : length (ticks)
                        labels_out{ii} = labels(ii, :);
                    end
                else
                    labels = repmat (labels, reps, 1); % repeat labels
                    for ii = 1 : length (ticks)
                        labels_out{ii} = labels{ii};
                    end
                end
            end
            
        end
    end

end
