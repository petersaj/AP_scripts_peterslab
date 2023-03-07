function runBonsaiWorkflow(workflowPath, addOptArray, bonsaiExePath, external)
%% addOptArray = {'field', 'val', 'field2', 'val2'};
%% check the input
switch nargin
    case 1
        bonsaiExePath = bonsaiPath(64);
        addOptArray = '';
        addFlag = 0;
        external = 1;
    case 2
        bonsaiExePath = bonsaiPath(64);
        addFlag = 1;
        external = 1;
    case 3
        addFlag = 1;
        external = 1;
    case 4
        if isempty(bonsaiExePath)
            bonsaiExePath = bonsaiPath(64);
        end
        addFlag = 1;
        
    otherwise
        error('Too many variables');
end



%% write the command
optSuffix = '-p:';
startFlag = '--start';
noEditorFlag = '--noeditor'; %use this instead of startFlag to start Bonsai without showing the editor
cmdsep = ' ';
if external
command = [['"' bonsaiExePath '"'] cmdsep ['"' workflowPath '"'] cmdsep startFlag];
else
command = [['"' bonsaiExePath '"']  cmdsep ['"' workflowPath '"'] cmdsep noEditorFlag];
end
ii = 1;
commandAppend = '';
if addFlag
    while ii < size(addOptArray,2)
        commandAppend = [commandAppend cmdsep optSuffix addOptArray{ii} '="' num2str(addOptArray{ii+1}) '"' ];
        ii = ii+2;
    end
end

if external
    command = [command commandAppend ' &'];
else
    command = [command commandAppend];
end
%% run the command
[sysecho] = system(command);
end


%% parse the addOptArray
function out = parseAddOptArray(addOptArray)
addOpt = cell(length(addOptArray)/2,2);
for ii =  1:length(addOpt)
    addOpt{ii,1} = addOptArray{2*ii-1}; %gets the propriety
    addOpt{ii,2} = addOptArray{2*ii}; %gets the value
end
out = addOpt;
end