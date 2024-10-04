function a = qa_extract_data_1chan(filename,sortcode,chan,sets)
% function a = qa_extract_data_1chan(filename,sortcode,chan)
%
% Extracts data from a NEV file and puts it in a simple struct with
% data and codes aligned to stim onset as well as saccade onset
%
% sortcode is optional list of valid sort codes (default is 1)
% chan is the electrode channel (default is 1)
%
%sets = {[165,315], [15:45:330]};

if (nargin < 2)
    sortcode = 1;
end

if (nargin < 3)
    chan = 1;
end

splitData = strsplit(filename, '_');
monkey = splitData{1};
session = splitData{2};
task_name = splitData{3};
chamber = splitData{4};

% if iscell(filename)
%     disp 'BROKEN CODE. NEEDS UPDATING'
%     pause
%     nev = readNEVMulti(filename);
% else
    [nev,~] = read_nev(filename); % Updated Feb 2024 to *not* use readNEV
% end

if (isempty(nev))
    error('NEV File Not Found');
end

% look for cases where channel is invalid
if isempty(find(nev(:,1) == chan))
    error(['*** No spikes on electrode ',num2str(chan)]);
end

% find digital codes and spikes
digvals = find(nev(:,1) == 0);
spkvals = find(nev(:,1) == chan);

% look for cases where there are no spikes with that sort code
if isempty(find(nev(spkvals,2) == sortcode))
    error(['*** No spikes with sort code ',num2str(sortcode),' on electrode ',num2str(chan)]);
end

% this is 1 electrode, so get rid of all codes other than 0
% (digital codes) or chan (spikes on channel 'chan')
nev = nev(sort([digvals;spkvals]),:);

% digital codes
dat = nev(nev(:,1) == 0,:);
dat = dat(dat(:,2) ~= 0,:);

tstarts = find(dat(:,2)==1);
tstops = find(dat(:,2)==255);
cnd=dat(tstarts+1,2)-2^15;
ncnd = length(unique(cnd));

%%%%%%%%%%%%%%%%% temporary add by Kendra to narrow down conditions
% Define the sets
[xs, ys] = ndgrid(sets{:});
cartProd = [xs(:) ys(:)];

% Define the number of groups and boundaries
num_groups = 8;
boundaries = linspace(0, 360, num_groups + 1);
angs = sets{2};

% Initialize cell array to store conditions
conditions = arrayfun(@(i) angs(angs >= boundaries(i) & angs < boundaries(i + 1)), 1:num_groups, 'UniformOutput', false);

% Initialize cell array to store groups
groups = arrayfun(@(i) find(ismember(cartProd(:,2), conditions{i})), 1:num_groups, 'UniformOutput', false);

% Update the cnd array based on the new groupings
new_cnd = arrayfun(@(x) find(cellfun(@(y) ismember(x, y), groups)), cnd);

cnd = new_cnd;
ncnd = length(unique(cnd));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EVENTS = cell X cnd X repeat
a = struct();

a.CND = unique(cnd);
a.RPTS = zeros(ncnd,1);

% in case the first '1' wasn't recorded
while tstops(1) < tstarts(1)
    tstops(1) = [];
end
% consider only completed trials
if length(tstarts) > length(tstops)
    tstarts(length(tstops)+1:end)=[];
end

for I=1:length(tstops)
    tcodes = dat(tstarts(I):tstops(I),:);
    % look for CORRECT code
    if (~isempty(find(tcodes(:,2)==5)))
        
        %rewtime = tcodes(find(tcodes(:,2)==73),3);
        % code 70 is TARG_ON
        stimtime=tcodes(find(tcodes(:,2)==70),3);
        stimtime = stimtime(1); % in case there are > 1 TARG_ON codes
        if isequal(task_name, 'dirmem')
            % code 141 is SACCADE (for reference, 3 is FIX_OFF)
            sactime=tcodes(find(tcodes(:,2)==141),3);
            sactime = sactime(1);
        else
            sactime = tcodes(find(tcodes(:,2)==70),3);
            sactime = sactime(1) + (110/1000);
        end
        % all nev codes
        tnev=nev(find(nev(:,3) >= tcodes(1,3) & nev(:,3) <= tcodes(end,3)),:);
        tnev(find(tnev(:,1)==0),:) = [];
        tnev((tnev(:,2) ~= sortcode),:) = [];
        % spike times
        st=tnev(:,3);
        % aligned on stim onset
        ttcodes = tcodes;
        ttcodes(:,3)=ttcodes(:,3)-stimtime;  
        a.STIMCODES{cnd(I),a.RPTS(cnd(I))+1} = ttcodes;
        a.STIMEVENTS{cnd(I),a.RPTS(cnd(I))+1} = st-stimtime;
        % aligned on sac onset (or fix offset)
        ttcodes = tcodes;
        ttcodes(:,3)=ttcodes(:,3)-sactime;  
        a.SACCODES{cnd(I),a.RPTS(cnd(I))+1} = ttcodes;
        a.SACEVENTS{cnd(I),a.RPTS(cnd(I))+1} = st-sactime;
        % increment repeats for that condition
        a.RPTS(cnd(I)) = a.RPTS(cnd(I))+1;
    end
end


x=1;

