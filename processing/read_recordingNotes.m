function [this_sess,filename] = read_recordingNotes(NOTES_PATH,EXPERIMENTER,MONKEY,SESSION)
    %dataFolder = '/Users/kendranoneman/OneDrive/DATA';
    %EXPERIMENTER = 'kendra';
    %MONKEY = 'scrappy';
    %SESSION = '0113a';
    
    recordings = readtable(NOTES_PATH);
    recordings.experimenter = categorical(recordings.experimenter);
    recordings.monkey = categorical(recordings.monkey);
    recordings.session_depth = categorical(recordings.session_depth);
    recordings.mapFile_name = categorical(recordings.mapFile_name );

    recordings.hemi = cellfun(@(q) eval(strrep(strrep(q, '“', '"'), '”', '"')), recordings.hemi, 'uni', 0);
    recordings.gridHole = cellfun(@(q) eval(strrep(strrep(q, '“', '"'), '”', '"')), recordings.gridHole, 'uni', 0);
    recordings.gtHeight_mm = cellfun(@(q) eval(strrep(strrep(q, '“', '"'), '”', '"')), recordings.gtHeight_mm, 'uni', 0);
    recordings.recordDepth_mm = cellfun(@(q) eval(strrep(strrep(q, '“', '"'), '”', '"')), recordings.recordDepth_mm, 'uni', 0);
    recordings.probeID = cellfun(@(q) eval(strrep(strrep(q, '“', '"'), '”', '"')), recordings.probeID, 'uni', 0);
    recordings.probeUse_num = cellfun(@(q) eval(strrep(strrep(q, '“', '"'), '”', '"')), recordings.probeUse_num, 'uni', 0);
    recordings.dead_ripChans = cellfun(@(q) eval(strrep(strrep(q, '“', '"'), '”', '"')), recordings.dead_ripChans, 'uni', 0);
    
    this_sess = recordings(recordings.experimenter==EXPERIMENTER & recordings.monkey==MONKEY & recordings.session_depth==SESSION,:);
    filename = sprintf('%s_%s_%s',EXPERIMENTER,MONKEY,SESSION);

end