%RSAC    Read SAC binary files.
%    RSAC('sacfile') reads in a SAC (seismic analysis code) binary
%    format file into a 3-column vector.
%    Column 1 contains time values.
%    Column 2 contains amplitude values.
%    Column 3 contains all SAC header information.
%    Default byte order is big-endian.  M-file can be set to default
%    little-endian byte order.
%
%    usage:  output = rsac('sacfile')
%
%    Examples:
%
%    KATH = rsac('KATH.R');
%    plot(KATH(:,1),KATH(:,2))
%
%    [SQRL, AAK] = rsac('SQRL.R','AAK.R');
%
%    by Michael Thorne (4/2004)   mthorne@asu.edu
%
%
% --------------------------------------------------------
% Nori
% 10/24/2013, add endian flag at variable, ='big-endian' or 'little-endian'
% data=rsac('big-endian',filename);
%
% --------------------------------------------------------
% Headers
%   1: DELTA (dt)
%   2: DEPMIN ?
%   3: DEPMAX ?
%   4: SCALE  ?
%   5: ODELTA ?
%   6: B (begin time from KZTIME)
%   7: E (end time from KZTIME)
%   8: O
%   9: A 
%  10: INTERNAL
%  11~20: T0~T9
%  21: F
%  22~31: RESP0~RESP9
%  32: STLA (station latitude)
%  33: STLO (station longitude)
%  34: STEL (station elevation (m))
%  35: STDP
%  36: EVLA
%  37: EVLO
%  38: EVEL
%  39: EVDP
%  40: MAG
%  41~50: USER0 ~ USER9
%  51: DIST
%  52: AZ
%  53: BAZ
%  58: CMPAZ (component azimuth)
%  59: CMPINC (component incident angle)
%  71: KZDATE (year)
%  72: KZDATE (julian day)
%  73: KZTIME (hour)
%  74: KZTIME (min)
%  75: KZTIME (sec)
%  76: KZTIME (sec after decimal point)
%  77: NVHDR (header version)
%  78: NORID
%  79: NEVID
%  80: NPTS (number of time sample)
%  86:
%  87:
% 106:
% 107:
% 108-306

%function [varargout] = rsac(varargin);
function [varargout] = rsac(endian,varargin);

for nrecs = 2:nargin
  nrecs=nrecs-1;

  sacfile = varargin{nrecs};

%for nrecs = 1:nargin

%  sacfile = varargin{nrecs};

%---------------------------------------------------------------------------
%    Default byte-order
%    endian  = 'big-endian' byte order (e.g., UNIX)
%            = 'little-endian' byte order (e.g., LINUX)

%endian = 'little-endian';
%endian = 'big-endian';

if strcmp(endian,'big-endian')
  fid = fopen(sacfile,'r','ieee-be'); 
elseif strcmp(endian,'little-endian')
  fid = fopen(sacfile,'r','ieee-le'); 
end

% read in single precision real header variables:
%---------------------------------------------------------------------------
for i=1:70
  h(i) = fread(fid,1,'single');
end

% read in single precision integer header variables:
%---------------------------------------------------------------------------
for i=71:105
  h(i) = fread(fid,1,'int32');
end


% Check header version = 6 and issue warning
%---------------------------------------------------------------------------
% If the header version is not NVHDR == 6 then the sacfile is likely of the
% opposite byte order.  This will give h(77) some ridiculously large
% number.  NVHDR can also be 4 or 5.  In this case it is an old SAC file
% and rsac cannot read this file in.  To correct, read the SAC file into
% the newest verson of SAC and w over.
% 
if (h(77) == 4 | h(77) == 5)
    message = strcat('NVHDR = 4 or 5. File: "',sacfile,'" may be from an old version of SAC.'); 
    error(message)
elseif h(77) ~= 6
    message = strcat('Current rsac byte order: "',endian,'". File: "',sacfile,'" may be of opposite byte-order.');
    error(message)
end

% read in logical header variables
%---------------------------------------------------------------------------
for i=106:110
  h(i) = fread(fid,1,'int32');
end

% read in character header variables
%---------------------------------------------------------------------------
for i=111:302
  h(i) = (fread(fid,1,'char'))';
end

% read in amplitudes
%---------------------------------------------------------------------------

YARRAY     = fread(fid,'single');
h(106)=1;
if h(106) == 1
  XARRAY = (linspace(h(6),h(7),h(80)))'; 
else
  error('LEVEN must = 1; SAC file not evenly spaced')
end 

% add header signature for testing files for SAC format
%---------------------------------------------------------------------------
h(303) = 77;
h(304) = 73;
h(305) = 75;
h(306) = 69;

% arrange output files
%---------------------------------------------------------------------------
OUTPUT(:,1) = XARRAY;
OUTPUT(:,2) = YARRAY;
OUTPUT(1:306,3) = h(1:306)';

%pad xarray and yarray with NaN if smaller than header field
if h(80) < 306
  OUTPUT((h(80)+1):306,1) = NaN;
  OUTPUT((h(80)+1):306,2) = NaN;
end

fclose(fid);

varargout{nrecs} = OUTPUT;

end
