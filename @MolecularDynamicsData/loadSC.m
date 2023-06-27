function loadSC(this, varargin)
    if nargin > 1
		filepath = varargin{1};
	else
		[dir, name, ext] = fileparts(this.StructureFile);
        filepath = fullfile(dir, "\sc.txt");
    end

	fid = fopen(filepath);
	tline = fgetl(fid);
% 	data = cell2table(cell(0,7), 'VariableNames', {'Chain1', 'Pos1', 'Res1', 'Chain2', 'Pos2', 'Res2', 'Frames'});
	data = cell(0,7);
	frame = 1;
	while ischar(tline)
		if strcmp(tline(1), '#') == 1
			[mat, tokens] = regexp(tline, '# (\d{1,8})', 'match', 'tokens');
			frame = str2double(tokens{1}{1}) + 1;
			disp(frame);
		else
			[mat, tokens] = regexp(tline, '(\w{1})(\d{1,3})(\w{1})-(\w{1})(\d{1,3})(\w{1})', 'match', 'tokens');
% 			ind = find(strcmp(data.Chain1, tokens{1}{1}) & data.Pos1 == str2double(tokens{1}{2}) & strcmp(data.Chain2, tokens{1}{4}) & data.Pos2 == str2double(tokens{1}{5}), 1, 'first');
			ind = [];
			pos1 = str2double(tokens{1}{2});
			pos2 = str2double(tokens{1}{5});
			if (pos1 >= 47 && pos1 <= 73) || (pos1 >= 177 && pos1 <= 203) || (pos2 >= 47 && pos2 <= 73) || (pos2 >= 177 && pos2 <= 203)
				if size(data, 1) > 0
					ind = find(find_str_cell(data(:,1), tokens{1}{1}) & cell2mat(data(:,2)) == pos1 & find_str_cell(data(:,4), tokens{1}{4}) & cell2mat(data(:,5)) == pos2, 1, 'first');
				end
				if isempty(ind)
					% New stab center
					data = [data; {tokens{1}{1}, str2double(tokens{1}{2}), tokens{1}{3}, tokens{1}{4}, str2double(tokens{1}{5}), tokens{1}{6}, frame}];
				else
					% Stab center already exists
					data{ind,7}(end+1) = frame;
				end
			end
		end

		tline = fgetl(fid);
	end
	fclose(fid);

	data = cell2table(data);
	data.Properties.VariableNames = {'Chain1', 'Pos1', 'Res1', 'Chain2', 'Pos2', 'Res2', 'Frames'};
	
	% Place residue with smaller number to first place
	for i = 1:height(data)
		if data.Pos2(i) < data.Pos1(i)
			chain1 = data.Chain1(i);
			pos1 = data.Pos1(i);
			res1 = data.Res1(i);
			data.Chain1(i) = data.Chain2(i);
			data.Pos1(i) = data.Pos2(i);
			data.Res1(i) = data.Res2(i);
			data.Chain2(i) = chain1;
			data.Pos2(i) = pos1;
			data.Res2(i) = res1;
		end
	end
	
	% Sort 
	data = sortrows(data, {'Pos1', 'Pos2', 'Chain1', 'Chain2'}, 'ascend');

	this.StabCenters = data;
end