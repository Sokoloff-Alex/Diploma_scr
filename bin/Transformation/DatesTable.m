function [] = DatesTable(varargin)
% function to print table of dates in different format
% 'yyyy-mm-dd     MJD        matlabTime  EpochNumber'
% to make date conversion for GMT faster
%
% Intut: StartDate, EndDate in [yyyy-mm-dd] format
% Outut: printed to stdout
%
% Alexandr Sokolov, KEG


if isempty(varargin)
    StartDate = datenum('2004-01-01');
    EndDate   = datenum('2017-01-01');
elseif length(varargin) == 2
    StartDate = datenum(varargin{1});
    EndDate   = datenum(varargin{2});
end

try
    MJD = 53005; % MJD of 2004-01-01
    MJD = MJD - 1;
    fprintf('yyyy-mm-dd     MJD        matlabTime  EpochNumber \n');
    EpochNr = 0;
    for t = StartDate:EndDate
        MJD = MJD + 1;
        EpochNr = EpochNr + 1;
        fprintf('%s     %s      %s      %s \n', datestr(t,'yyyy-mm-dd'), num2str(MJD), num2str(t), num2str(EpochNr));   
    end
catch
        disp('use correct Start/End dates format')
end

end