classdef MolecularDynamicsData < handle
    properties (SetObservable)
        % Structure file (.pdb or .cms)
        StructureFile string
        % Trajectory file (.dcd or .dtr)
        TrajectoryFile string
		
		% Frame interval in ps
		FrameInterval
		
		% Sequence and region assignation
		Sequence table
		
		% VMD TCL script text
        VMDscript struct
    end
    properties (Transient, SetObservable)
		% Stabilization centers
		StabCenters table
 
		% Position of C-alpha of selected residues
		ResiduePositionCA table
		ResiduePositionS table
		
		% Torsion angles of proline residues
		TorsionAngles table

		% RMSD compared to a reference
		RMSD table
		
		% Center of mass of hemichannels
		CenterOfMass

		% Distance between different atoms
		Distances struct
		DistancesCA table
		DistancesS table
		DistancesSacceptor table
		DistancesSdonor table
		
		% Hydrogen bonds
		Hbonds table

		% pKa values
		pKa table
    end
    properties (SetAccess = private)
		%
    end
    properties (Transient, SetAccess = private, SetObservable)
        % 
    end
    properties (Dependent, SetObservable)
		% Folder containing the data
		DataFolder string

		% Number of frames
		FrameNum double
		
		% Frame times (in ns)
        Time double
		
		% Filtered stab centers
		StabCentersFrequent
		StabCentersExtracell
		StabCentersExtracellFrequent
		StabCentersDisulfide
		StabCentersTransjunctional
    end
    
    methods (Static)
        %
	end
   
    methods
        % Set methods
     
        % Get methods
		function value = get.DataFolder(this)
			[path, name, ext] = fileparts(this.StructureFile);
			value = path;
		end
		
		function value = get.FrameNum(this)
			if ~isempty(this.StabCenters) && ~isempty(this.Distances)
				value = double(min([max([this.StabCenters.Frames{:}]) max(this.Distances.SS.Frame)]));
			elseif ~isempty(this.StabCenters)
				value = double(max([this.StabCenters.Frames{:}]));
			elseif ~isempty(this.Distances)
				value = double(max(this.Distances.SS.Frame));
			elseif ~isempty(this.ResiduePositionCA)
				value = double(max(this.ResiduePositionCA.frame)+1);
			end
		end
        
		function value = get.Time(this)
			value = (1:this.FrameNum)' * this.FrameInterval / 1000;
		end
        
		function value = get.StabCentersFrequent(this)
			occurence = cellfun(@length, this.StabCenters.Frames);
			filter = occurence > this.FrameNum * 0.02;
			value = this.StabCenters(filter, :);
		end
        
		function value = get.StabCentersExtracell(this)
			ec = this.Sequence.Region(this.Sequence.Pos) < 0;
			filter = ec(this.StabCenters.Pos1) | ec(this.StabCenters.Pos2);
			value = this.StabCenters(filter, :);
		end
        
		function value = get.StabCentersExtracellFrequent(this)
			ec = this.Sequence.Region(this.Sequence.Pos) < 0;
			occurence = cellfun(@length, this.StabCenters.Frames);
			filter = (ec(this.StabCenters.Pos1) | ec(this.StabCenters.Pos2)) & occurence > this.FrameNum * 0.02;
			value = this.StabCenters(filter, :);
		end
        
		function value = get.StabCentersDisulfide(this)
			filter = this.StabCenters.Res1 == "C" & this.StabCenters.Res2 == "C";
			value = this.StabCenters(filter, :);
		end
        
		function value = get.StabCentersTransjunctional(this)
			hc1 = {'A', 'B', 'C', 'D', 'E', 'F'};
			hc2 = {'G', 'H', 'I', 'J', 'K', 'L'};
			filter = (cellfun(@(x) sum(find_str_cell(hc1, x)) > 0, this.StabCentersExtracellFrequent.Chain1) & cellfun(@(x) sum(find_str_cell(hc2, x)) > 0, this.StabCentersExtracellFrequent.Chain2)) | (cellfun(@(x) sum(find_str_cell(hc2, x)) > 0, this.StabCentersExtracellFrequent.Chain1) & cellfun(@(x) sum(find_str_cell(hc1, x)) > 0, this.StabCentersExtracellFrequent.Chain2));
			value = this.StabCentersFrequent(filter, :);
		end
        
		function value = get.VMDscript(this)
			% Load structure and trajectory files
% 			script.load = "mol new {" + this.StructureFile + "} type {mae} first 0 last -1 step 1 waitfor 1";
% 			script.load(end+1,1) = "mol addfile {" +  this.TrajectoryFile + "} type {dtr} first 0 last -1 step 100 waitfor 1 [molinfo top]";
% 			script.load(end+1,1) = "animate delete beg 0 end 0 skip 0 [molinfo top]";
			
% 			script.createMaterials = strings(height(this.StabCentersExtracellFrequent)*2, 1);
% 			for i = 1:height(this.StabCentersExtracellFrequent)
% 				aa = strcat(this.StabCentersExtracellFrequent.Chain1{i}, num2str(this.StabCentersExtracellFrequent.Pos1(i)));
% 				script.createMaterials(i*2-1,1) = strcat("material add Transparent", aa, " copy Transparent");
% 				aa = strcat(this.StabCentersExtracellFrequent.Chain2{i}, num2str(this.StabCentersExtracellFrequent.Pos2(i)));
% 				script.createMaterials(i*2,1) = strcat("material add Transparent", aa, " copy Transparent");
% 			end

			% Reset			
			script.sc = strings(0);
			script.sc{end+1,1} = 'set sel [atomselect top "protein"]';
			script.sc{end+1,1} = '$sel set beta 0;';
			script.sc{end+1,1} = 'for {set iFrame 0} {$iFrame < [molinfo top get numframes]} {incr iFrame} {';
			script.sc{end+1,1} = '$sel frame $iFrame; $sel set user 0;';
			script.sc{end+1,1} = '}';
			script.sc{end+1,1} = '$sel delete';
				
			vmdframereduction = 10;
			stabCenters = this.StabCentersExtracell;
			sc = unique([strcat(stabCenters.Chain1, num2str(stabCenters.Pos1)); strcat(stabCenters.Chain2, num2str(stabCenters.Pos2))]);
			for iSC = 1:length(sc)
				aa = sc{iSC}(1);
				pos = str2double(sc{iSC}(2:end));
				f = stabCenters.Frames((find_str_cell(stabCenters.Chain1, aa) == 1 & stabCenters.Pos1 == pos) | (find_str_cell(stabCenters.Chain2, aa) == 1 & stabCenters.Pos2 == pos))';
				f = unique(cell2mat(f));
				if ~isempty(f)
					script.sc{end+1,1} = sprintf('set sel [atomselect top "chain %s and resid %d"]', aa, pos);
					% Relative occurrence of this SC during the whooe MD
					script.sc{end+1,1} = sprintf('$sel set beta %d;', length(f) / this.FrameNum * 100);
					f = unique(round(f ./ vmdframereduction));
					for iFrame = f
						script.sc{end+1,1} = sprintf('$sel frame %d; $sel set user 100', iFrame);
					end
					script.sc{end+1,1} = '$sel delete';
				end
			end
			
			value = script;
		end
        
        % Simple methods
		function VMDwrite(this)
% 			fprintf(fid, '%s\n', this.VMDscript.load);
% 			for i = 1:1000:length(this.VMDscript.sc)
				fid = fopen("D:\MD\matlab.tcl", "wt");
				fprintf(fid, '%s\n', this.VMDscript.sc);
				fclose(fid);
% 			end
		end
		
		function VMDrun(this)
			system(strcat(char({'"c:\Program Files (x86)\University of Illinois\VMD\vmd.exe"'}), char({' "play "D:\MD\matlab.tcl"'})));
		end
		
		function loadCenterOfMass(this)
			filepath = fullfile(this.DataFolder, "centerofmass.txt");
			fid = fopen(filepath);
			ts = textscan(fid,'%s %u %u %f %f %f');
			fclose(fid);
			
			this.CenterOfMass = table(string(ts{1}), ts{3}, ts{4}, ts{5}, ts{6});
			this.CenterOfMass.Properties.VariableNames = {'HC' 'frame' 'x' 'y' 'z'};
		end
		
		function loadDistancesSacceptor(this)
			filepath = fullfile(this.DataFolder, "distSacceptor.txt");
			fid = fopen(filepath);
			ts = textscan(fid,'%u %s %s %u %s %s %s %u %s %f');
			fclose(fid);
			
			this.DistancesSacceptor = table(ts{1}, string(ts{2}), string(ts{3}), ts{4}, string(ts{5}), string(ts{6}), string(ts{7}), ts{8}, string(ts{9}), ts{10});
			this.DistancesSacceptor.Properties.VariableNames = {'Frame' 'Chain1' 'Res1' 'Pos1' 'Atom1' 'Chain2' 'Res2' 'Pos2' 'Atom2' 'Dist'};
		end
		
		function loadDistancesSdonor(this)
			filepath = fullfile(this.DataFolder, "distSdonor.txt");
			fid = fopen(filepath);
			ts = textscan(fid,'%u %s %s %u %s %s %s %u %s %f');
			fclose(fid);
			
			this.DistancesSdonor = table(ts{1}, string(ts{2}), string(ts{3}), ts{4}, string(ts{5}), string(ts{6}), string(ts{7}), ts{8}, string(ts{9}), ts{10});
			this.DistancesSdonor.Properties.VariableNames = {'Frame' 'Chain1' 'Res1' 'Pos1' 'Atom1' 'Chain2' 'Res2' 'Pos2' 'Atom2' 'Dist'};
		end
		
		function loadHbondsTransGJ(this)
			filepath = fullfile(this.DataFolder, "transGJHbond, heavy-H.txt");
			fid = fopen(filepath);
			ts = textscan(fid,'%u %s %s %u %s %s %s %u %s %f');
			fclose(fid);
			
			this.Hbonds = table(ts{1}, string(ts{2}), string(ts{3}), ts{4}, string(ts{5}), string(ts{6}), string(ts{7}), ts{8}, string(ts{9}), ts{10});
			this.Hbonds.Properties.VariableNames = {'Frame' 'Chain1' 'Res1' 'Pos1' 'Atom1' 'Chain2' 'Res2' 'Pos2' 'Atom2' 'Dist'};
		end
		
		function loadHbondsCYS(this)
			filepath = fullfile(this.DataFolder, "CYS-55-58Hbond.txt");
			fid = fopen(filepath);
			ts = textscan(fid,'%u %s %s %u %s %s %s %u %s %f');
			fclose(fid);
			
			this.Hbonds = table(ts{1}, string(ts{2}), string(ts{3}), ts{4}, string(ts{5}), string(ts{6}), string(ts{7}), ts{8}, string(ts{9}), ts{10});
			this.Hbonds.Properties.VariableNames = {'Frame' 'Chain1' 'Res1' 'Pos1' 'Atom1' 'Chain2' 'Res2' 'Pos2' 'Atom2' 'Dist'};
		end
		
		function loadResPos(this)
			% Positions of CA atoms
			filepath = fullfile(this.DataFolder, "resposCA.txt");
			fid = fopen(filepath);
			ts = textscan(fid,'%s %u %u {%f %f %f}');
			fclose(fid);
			
			this.ResiduePositionCA = table(string(ts{1}), ts{2}, ts{3}, ts{4}, ts{5}, ts{6});
			this.ResiduePositionCA.Properties.VariableNames = {'chain' 'residue' 'frame' 'x' 'y' 'z'};
			
			% Positions of S atoms
			filepath = fullfile(this.DataFolder, "resposS.txt");
  			fid = fopen(filepath);
			ts = textscan(fid,'%s %u %u {%f %f %f}');
			fclose(fid);
			
			this.ResiduePositionS = table(string(ts{1}), ts{2}, ts{3}, ts{4}, ts{5}, ts{6});
			this.ResiduePositionS.Properties.VariableNames = {'chain' 'residue' 'frame' 'x' 'y' 'z'};
			this.calculateDistancesS;
		end
		
		function loadRMSD(this)
			filepath = fullfile(this.DataFolder, "rmsd.txt");
			fid = fopen(filepath);
			ts = textscan(fid,'%u %f %f');
			fclose(fid);
			
			this.RMSD = table(ts{1}, ts{2}, ts{3});
			this.RMSD.Properties.VariableNames = {'frame' 'rmsd2start' 'rmsd2end'};
		end
		
		function loadRMSDres(this)
			filepath = fullfile(this.DataFolder, "RMSDres.txt");
			fid = fopen(filepath);
			ts = textscan(fid,'%s %u %u %f %f');
			fclose(fid);
			
			this.RMSD = table(string(ts{1}), ts{2}, ts{3}, ts{4}, ts{5});
			this.RMSD.Properties.VariableNames = {'chain' 'residue' 'frame' 'rmsd2start' 'rmsd2end'};
		end
		
		function loadpKa(this)
			files = dir(strcat(this.DataFolder, "/propka"));
			text = fileread(fullfile(files(3).folder, files(3).name));
			p = struct2table(regexp(text, '   (?<Res>[A-Z]{3}) {1,3}(?<Pos>\d{1,3}) (?<Chain>[A-Z]) +(?<pKa>[0-9.-]+)', 'names'));
			p.Pos = str2double(p.Pos);
			p.pKa = str2double(p.pKa);
			for iFile = 4:size(files, 1)
				disp(files(iFile).name);
				text = fileread(fullfile(files(iFile).folder, files(iFile).name));
				temp = str2double(struct2table(regexp(text, '   [A-Z]{3} {1,3}\d{1,3} [A-Z] +(?<pKa>[0-9.-]+)', 'names')').pKa);
				p.pKa = horzcat(p.pKa, temp);
			end
			
			this.pKa = p;
% 			this.DistancesSdonor.Properties.VariableNames = {'Frame' 'Chain1' 'Res1' 'Pos1' 'Atom1' 'Chain2' 'Res2' 'Pos2' 'Atom2' 'Dist'};
		end
		
		function loadTorsionAngles(this)
			filepath = fullfile(this.DataFolder, "torsionangles.txt");
			fid = fopen(filepath);
			ts = textscan(fid,'%s %u %u %f %f %f %f');
			fclose(fid);
			
			this.TorsionAngles = table(string(ts{1}), ts{2}, ts{3}, ts{4}, ts{5}, ts{6}, ts{7});
			this.TorsionAngles.Properties.VariableNames = {'chain' 'residue' 'frame' 'omegaprev' 'phi' 'psi' 'omega'};
		end
		
		function plotCenterOfMass(this)
			figure;
			plot(this.CenterOfMass.frame(this.CenterOfMass.HC == "HC1"), this.CenterOfMass.z(this.CenterOfMass.HC == "HC1"));
			hold on;
			plot(this.CenterOfMass.frame(this.CenterOfMass.HC == "HC2"), this.CenterOfMass.z(this.CenterOfMass.HC == "HC2"));
			plot(this.CenterOfMass.frame(this.CenterOfMass.HC == "HC2"), abs(this.CenterOfMass.z(this.CenterOfMass.HC == "HC1") - this.CenterOfMass.z(this.CenterOfMass.HC == "HC2")));
		end
		
		function plotHbonds(this)
% 			hbonds = this.Hbonds(this.Hbonds.Dist <= 2,:);
			hbonds = this.Hbonds(this.Hbonds.Frame > 0,:);
			hbondsall = [];
			hbondsallTickLabel = string.empty;
			hbondsatoms = [];
			hbondsatomsTickLabel = string.empty;
			uPos1 = unique(hbonds.Pos1);
			uPos2 = unique(hbonds.Pos2);
			for i = 1:length(uPos1)
				iPos1 = uPos1(i);
% 				uPos2 = uPos2(uPos2 >= iPos1);
				for j = 1:length(uPos2)
					iPos2 = uPos2(j);
					hbondsPair = hbonds((hbonds.Pos1 == iPos1 & hbonds.Pos2 == iPos2) | (hbonds.Pos1 == iPos2 & hbonds.Pos2 == iPos1),:);
					if height(hbondsPair) > 0
						% Individual pairs on individual chains
						hbondsPair = sortrows(hbondsPair, [4 8 5 9 2 1]);
						[u,ia,iu] = unique(hbondsPair(:, 2:9), 'rows', 'stable');
						resovertime = false(height(u), 10000);

						for iPair = 1:height(u)
							resovertime(iPair, hbondsPair.Frame(iu == iPair)) = true;
						end
% 
						% Individual pairs on individual chains
						[u,ia,iu] = unique(u(:, [2 3 4 6 7 8]), 'rows', 'stable');

						for iPair = 1:height(u)
							hbondsatoms(end+1,:) = sum(resovertime(iu == iPair,1:10000), 1);
							hbondsatomsTickLabel(size(hbondsatoms, 1), 1) = strcat(u.Res1(iPair), string(u.Pos1(iPair)), ' "', u.Atom1(iPair), '"', " - ", u.Res2(iPair), string(u.Pos2(iPair)), ' "', u.Atom2(iPair), '"');
						end

						% All pairs on all chains
						hbondsall(end+1,:) = sum(resovertime(:,1:10000), 1);
						hbondsallTickLabel(size(hbondsall, 1), 1) = strcat(string(iPos1), aa3to1(u.Res1(1)), " - ", string(iPos2), aa3to1(u.Res2(1)));
					end
				end
			end
			
			% Individual pairs on individual chains
			figure('Name', inputname(1));
			imagesc(hbondsatoms);
			c = colormap(jet);
			c(1,3) = 0;
			colormap(c);
			set(gca, 'CLim', [0 12]);
			xlabel('Time (ns)');
			set(gca, 'YTick', 1:length(hbondsatomsTickLabel));
			set(gca, 'YTickLabel', hbondsatomsTickLabel);
			spaceplots;
			colorbar;
			
			% All pairs on all chains
			figure('Name', inputname(1));
			imagesc(hbondsall);
			c = colormap(jet);
			c(1,3) = 0;
			colormap(c);
			ax = gca;
			ax.CLim = [0 12];
			xlabel('Time (ns)');
			set(gca, 'YTickLabel', hbondsallTickLabel);
			spaceplots;
			colorbar;
		end
		
		function plotStabCenters(this, stabCenters)
			resovertime = false(height(stabCenters), this.FrameNum);
			for iRes = 1:height(stabCenters)
				resovertime(iRes, stabCenters.Frames{iRes}) = true;
			end
			
			transjunction = (ismember(stabCenters.Chain1, {'A', 'B', 'C', 'D', 'E', 'F'}) & ismember(stabCenters.Chain2, {'G', 'H', 'I', 'J', 'K', 'L'})) | (ismember(stabCenters.Chain2, {'A', 'B', 'C', 'D', 'E', 'F'}) & ismember(stabCenters.Chain1, {'G', 'H', 'I', 'J', 'K', 'L'}));
			transjunctionNum = sum(transjunction);
			intersubunit = ~transjunction & strcmp(stabCenters.Chain1, stabCenters.Chain2) == 0;
			intersubunitNum = sum(intersubunit);
			interloop = strcmp(stabCenters.Chain1, stabCenters.Chain2) == 1 & abs(stabCenters.Pos1 - stabCenters.Pos2) > 50;
			interloopNum = sum(interloop);
			
			if transjunctionNum > 0
				figure('Name', inputname(1));
				imagesc(this.Time, 1:height(stabCenters(transjunction,:)), resovertime(transjunction,:));
				title(strcat(inputname(1), ' stab. centers, transjunction'), 'Interpreter', 'none');
				cmapb = ones(2,3);
				cmapb(2,:) = [0 0 1];
				colormap(cmapb);
				xlabel('Time (ns)');
				set(gca, 'Ytick', 1:height(stabCenters(transjunction,:)), 'YTickLabel', strcat(stabCenters.Chain1(transjunction), num2str(stabCenters.Pos1(transjunction)), stabCenters.Res1(transjunction), '-', stabCenters.Chain2(transjunction), num2str(stabCenters.Pos2(transjunction)), stabCenters.Res2(transjunction)));
			end
			
			if intersubunitNum > 0
				figure('Name', inputname(1));
				imagesc(this.Time, 1:height(stabCenters(intersubunit,:)), resovertime(intersubunit,:));
				title(strcat(inputname(1), ' stab. centers, intersubunit'), 'Interpreter', 'none');
				cmapb = ones(2,3);
				cmapb(2,:) = [0 1 0];
				colormap(cmapb);
				xlabel('Time (ns)');
				set(gca, 'Ytick', 1:height(stabCenters(intersubunit,:)), 'YTickLabel', strcat(stabCenters.Chain1(intersubunit), num2str(stabCenters.Pos1(intersubunit)), stabCenters.Res1(intersubunit), '-', stabCenters.Chain2(intersubunit), num2str(stabCenters.Pos2(intersubunit)), stabCenters.Res2(intersubunit)));
			end
			
			if interloopNum > 0
				sc = stabCenters(interloop,:);
				rov = resovertime(interloop,:);
				for i = 1:50:height(sc)
					figure('Name', inputname(1));
					sci = sc(i:min(i+49,height(sc)),:);
					rovi = rov(i:min(i+49,height(sc)),:);
					imagesc(this.Time, 1:height(sci), rovi);
					title(strcat(inputname(1), ' stab. centers, interloop'), 'Interpreter', 'none');
					cmapr = ones(2,3);
					cmapr(2,:) = [1 0 0];
					colormap(cmapr);
					xlabel('Time (ns)');
					set(gca, 'Ytick', 1:height(sci), 'YTickLabel', strcat(sci.Chain1, num2str(sci.Pos1), sci.Res1, '-', sci.Chain2, num2str(sci.Pos2), sci.Res2));
				end
			end
		end
		
		function plotStabCentersRes(this, stabCenters)
			stabCentersRes = table;
			uPos1 = unique(stabCenters.Pos1);
			for iPos1 = 1:length(uPos1)
				iRes1 = unique(stabCenters.Res1(stabCenters.Pos1 == uPos1(iPos1)));
				uPos2 = unique(stabCenters.Pos2(stabCenters.Pos1 == uPos1(iPos1)));
				for iPos2 = 1:length(uPos2)
					iRes2 = unique(stabCenters.Res2(stabCenters.Pos2 == uPos2(iPos2)));
					sc = stabCenters(stabCenters.Pos1 == uPos1(iPos1) & stabCenters.Pos2 == uPos2(iPos2), :);
					f = zeros(1, this.FrameNum);
					for i = 1:height(sc)
						ft = zeros(1, this.FrameNum);
						ft(1, sc.Frames{i}) = 1;
						f = f + ft;
					end
					t = table(uPos1(iPos1), string(iRes1{1}), uPos2(iPos2), string(iRes2{1}), f, 'VariableNames', {'Pos1', 'Res1', 'Pos2', 'Res2', 'Frames'});
					stabCentersRes = [stabCentersRes; t];
				end
			end
			
			resovertime = false(height(stabCentersRes), this.FrameNum);
			for iRes = 1:height(stabCentersRes)
				resovertime(iRes, stabCentersRes.Frames{iRes}) = true;
			end
			
			transjunction = (ismember(stabCentersRes.Chain1, {'A', 'B', 'C', 'D', 'E', 'F'}) & ismember(stabCentersRes.Chain2, {'G', 'H', 'I', 'J', 'K', 'L'})) | (ismember(stabCentersRes.Chain2, {'A', 'B', 'C', 'D', 'E', 'F'}) & ismember(stabCentersRes.Chain1, {'G', 'H', 'I', 'J', 'K', 'L'}));
			transjunctionNum = sum(transjunction);
			intersubunit = ~transjunction & strcmp(stabCentersRes.Chain1, stabCentersRes.Chain2) == 0;
			intersubunitNum = sum(intersubunit);
			interloop = strcmp(stabCentersRes.Chain1, stabCentersRes.Chain2) == 1 & abs(stabCentersRes.Pos1 - stabCentersRes.Pos2) > 50;
			interloopNum = sum(interloop);
			
			if transjunctionNum > 0
				figure('Name', inputname(1));
				imagesc(this.Time, 1:height(stabCentersRes(transjunction,:)), resovertime(transjunction,:));
				title(strcat(inputname(1), ' stab. centers, transjunction'));
				cmapb = ones(2,3);
				cmapb(2,:) = [0 0 1];
				colormap(cmapb);
				xlabel('Time (ns)');
				set(gca, 'Ytick', 1:height(stabCentersRes(transjunction,:)), 'YTickLabel', strcat(stabCentersRes.Chain1(transjunction), num2str(stabCentersRes.Pos1(transjunction)), stabCentersRes.Res1(transjunction), '-', stabCentersRes.Chain2(transjunction), num2str(stabCentersRes.Pos2(transjunction)), stabCentersRes.Res2(transjunction)));
			end
			
			if intersubunitNum > 0
				figure('Name', inputname(1));
				imagesc(this.Time, 1:height(stabCentersRes(intersubunit,:)), resovertime(intersubunit,:));
				title(strcat(inputname(1), ' stab. centers, intersubunit'));
				cmapb = ones(2,3);
				cmapb(2,:) = [0 1 0];
				colormap(cmapb);
				xlabel('Time (ns)');
				set(gca, 'Ytick', 1:height(stabCentersRes(intersubunit,:)), 'YTickLabel', strcat(stabCentersRes.Chain1(intersubunit), num2str(stabCentersRes.Pos1(intersubunit)), stabCentersRes.Res1(intersubunit), '-', stabCentersRes.Chain2(intersubunit), num2str(stabCentersRes.Pos2(intersubunit)), stabCentersRes.Res2(intersubunit)));
			end
			
			if interloopNum > 0
				sc = stabCentersRes(interloop,:);
				rov = resovertime(interloop,:);
				for i = 1:50:height(stabCentersRes(interloop,:))
					figure('Name', inputname(1));
					sci = sc(i:min(i+49,height(sc)),:);
					rovi = rov(i:min(i+49,height(sc)),:);
					imagesc(this.Time, 1:height(sci), rovi);
					title(strcat(inputname(1), ' stab. centers, interloop'));
					cmapr = ones(2,3);
					cmapr(2,:) = [1 0 0];
					colormap(cmapr);
					xlabel('Time (ns)');
					set(gca, 'Ytick', 1:height(sci), 'YTickLabel', strcat(sci.Chain1, num2str(sci.Pos1), sci.Res1, '-', sci.Chain2, num2str(sci.Pos2), sci.Res2));
				end
			end
		end
		
		function plotDist(this)
			res = string(unique([this.Distances.SS.Res1; this.Distances.SS.Res2]));
			for iRes = 1%:length(res)
				figure;
				plot(this.Time, this.Distances.SSminValue(iRes,:));
				title(res(iRes));
				xlabel('Time (ns)');
				ylabel('Distance (A)');
				partners = unique(this.Distances.SSminPartner(iRes,:));
				partners(isnan(partners)) = [];
				ylim([3 6 + 0.2 * (length(partners) + 1)]);
				hold on;
				count = 1;
				for iPartner = partners
					y = 6 + 0.2 * count * (this.Distances.SSminPartner(iRes,:) == iPartner);
					y(y == 6) = nan;
					scatter(this.Time, y, 'r.');
					text(101, 6 + 0.2 * count, res(iPartner));
					count = count + 1;
				end
			end
		end
        
		function plotDistS(this)
			res = string(unique([this.Distances.SS.Res1; this.Distances.SS.Res2]));
			for iRes = 1%:length(res)
				figure;
				plot(this.Time, this.Distances.SSminValue(iRes,:));
				title(res(iRes));
				xlabel('Time (ns)');
				ylabel('Distance (A)');
				partners = unique(this.Distances.SSminPartner(iRes,:));
				partners(isnan(partners)) = [];
				ylim([3 6 + 0.2 * (length(partners) + 1)]);
				hold on;
				count = 1;
				for iPartner = partners
					y = 6 + 0.2 * count * (this.Distances.SSminPartner(iRes,:) == iPartner);
					y(y == 6) = nan;
					scatter(this.Time, y, 'r.');
					text(101, 6 + 0.2 * count, res(iPartner));
					count = count + 1;
				end
			end
		end
        
		function plotDistSacceptor(this)
			pos = unique([this.DistancesSacceptor.Pos1 this.DistancesSacceptor.Pos2], 'row');
			for iPos = 1:size(pos, 1)
				dist = this.DistancesSacceptor(this.DistancesSacceptor.Pos1 == pos(iPos, 1) & this.DistancesSacceptor.Pos2 == pos(iPos, 2) & this.DistancesSacceptor.Dist <= 3.2, :);
				if ~isempty(dist)
					figure;
					chains = unique(dist.Chain1);
					l = {};
					for iChain = 1:size(chains, 1)
						distChain = dist(dist.Chain1 == chains(iChain),:);
						scatter(distChain.Frame * this.FrameInterval / 1000, distChain.Dist, 'o');
						hold on;
						l{end+1,1} = strcat("Subunit ", chains(iChain));
					end
					title(strcat(dist.Res1(1), string(dist.Pos1(1)), "-", dist.Res2(1), string(dist.Pos2(1))));
					xlim([0 100]);
					ylim([2.5 3.2]);
					xlabel('Time (ns)');
					ylabel('Distance (A)');
					legend(l);
				end
			end
		end
        
		function plotDistSdonor(this)
			pos = unique([this.DistancesSdonor.Pos1 this.DistancesSdonor.Pos2], 'row');
			for iPos = 1:size(pos, 1)
				dist = this.DistancesSdonor(this.DistancesSdonor.Pos1 == pos(iPos, 1) & this.DistancesSdonor.Pos2 == pos(iPos, 2) & this.DistancesSdonor.Dist <= 3.2, :);
				if ~isempty(dist)
					figure;
					chains = unique(dist.Chain1);
					l = {};
					for iChain = 1:size(chains, 1)
						distChain = dist(dist.Chain1 == chains(iChain),:);
						scatter(distChain.Frame * this.FrameInterval / 1000, distChain.Dist, 'o');
						hold on;
						l{end+1,1} = strcat("Subunit ", chains(iChain));
					end
					title(strcat(dist.Res1(1), string(dist.Pos1(1)), "-", dist.Res2(1), string(dist.Pos2(1))));
					xlim([0 100]);
					ylim([2.5 3.2]);
					xlabel('Time (ns)');
					ylabel('Distance (A)');
					legend(l);
				end
			end
		end
        
		function plotRMSD(this)
			chain = unique(this.RMSD.chain);
			for iChain = 1:length(chain)
				data2start = this.RMSD.rmsd2start(this.RMSD.chain == chain(iChain));
				data2end = this.RMSD.rmsd2end(this.RMSD.chain == chain(iChain));
				if chain(iChain) == "EL"
					figure;
					plot((1:max(double(this.RMSD.frame))+1) / this.FrameInterval, reshape(data2start, [], max(this.RMSD.frame)+1));
					hold on;
					plot((1:max(double(this.RMSD.frame))+1) / this.FrameInterval, reshape(data2end, [], max(this.RMSD.frame)+1));
					title(this.DataFolder);
					xlabel("Time (ns)");
					ylabel("RMSD");
				else
					figure;
					imagesc(reshape(data2start, [], max(this.RMSD.frame)+1));
					title(strcat("chain ", chain(iChain)));
					xlabel("Frame");
					ylabel("Residue");
				end
			end
		end
         
		function plotRMSDres(this)
			res = unique(this.RMSD.residue);
			chain = unique(this.RMSD.chain);
			for iRes = 1:length(res)
				figure;
				for iChain = 1:length(chain)
					if chain(iChain) ~= "EL"
						plot((1:max(double(this.RMSD.frame))+1) / this.FrameInterval, this.RMSD.rmsd2start(this.RMSD.chain == chain(iChain) & this.RMSD.residue == res(iRes)));
						hold on;
						plot((1:max(double(this.RMSD.frame))+1) / this.FrameInterval, this.RMSD.rmsd2end(this.RMSD.chain == chain(iChain) & this.RMSD.residue == res(iRes)));
						xlabel("Time (ns)");
						ylabel("RMSD");
					end
				end
				title(strcat("residue ", num2str(res(iRes))));
				figure;
				data2start = this.RMSD.rmsd2start(this.RMSD.residue == res(iRes));
				plot((1:max(double(this.RMSD.frame))+1) / this.FrameInterval, mean(reshape(data2start, [], max(this.RMSD.frame)+1)));
			end
		end
         
		function plotRMSDstart2end(this)
			res = unique(this.RMSD.residue(strlength(this.RMSD.chain) == 1));
			chain = unique(this.RMSD.chain(strlength(this.RMSD.chain) == 1));
			for iChain = 1:length(chain)
				figure;
				for iRes = 1:length(res)
					plot(res, this.RMSD.rmsd(this.RMSD.chain == chain(iChain) & this.RMSD.frame == 0), 'ro');
					hold on;
					plot(res, this.RMSD.rmsd(this.RMSD.chain == chain(iChain) & this.RMSD.frame == max(this.RMSD.frame)), 'bo');
					xlabel("Residue");
					ylabel("RMSD");
				end
				title(strcat("chain ", chain(iChain)));
			end
		end
		
        % Complex methods
		result = calculateDistances(this, resids)
		calculateDistancesCA(this)
		calculateDistancesS(this)
		result = findNearest(this, resid)
        loadSC(this, filepath)
		reportAngles(this)
		reportAngleSCComparison(this)
		reportDistances(this)
		reportSC(this)
    end
end