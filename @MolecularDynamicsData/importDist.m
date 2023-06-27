function importDist(this, varargin)
    if nargin > 1
		filepath = varargin{1};
	else
		[dir, name, ext] = fileparts(this.StructureFile);
        filepath = fullfile(dir, "\distSS.txt");
	end
	
	% Load data
	data = cell2table(readcell(filepath));
	data.Properties.VariableNames = {'Frame', 'Res1', 'Res2', 'Dist'};

	this.Distances(1).SS = data;

	% Calculate minimum distance and determine disulfide partner of each Cys
	res = string(unique([this.Distances.SS.Res1; this.Distances.SS.Res2]));
	mindist = nan(length(res), this.FrameNum);
	partner = nan(length(res), this.FrameNum);
	for iRes = 1:length(res)
		t = this.Distances.SS(this.Distances.SS.Res1 == res(iRes) | this.Distances.SS.Res2 == res(iRes),:);
		for iFrame = 1:this.FrameNum
			if sum(t.Frame == iFrame) > 0
				mindist(iRes, iFrame) = min(t.Dist(t.Frame == iFrame));
				p = [t.Res1(t.Dist == mindist(iRes, iFrame) & t.Frame == iFrame) t.Res2(t.Dist == mindist(iRes, iFrame) & t.Frame == iFrame)];
				p = setdiff(p, res(iRes));
				partner(iRes, iFrame) = find(strcmp(res, p), 1);
			end
		end
	end
	
	this.Distances(1).SSminValue = mindist;
	this.Distances(1).SSminPartner = partner;
end