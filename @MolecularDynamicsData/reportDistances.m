function reportDistances(this)
	import mlreportgen.report.* 
	import mlreportgen.dom.* 
	
	rpt = Report(strcat("Distances, ", inputname(1)), "docx");

	add(rpt, Heading1("Összefoglalás"));
	add(rpt, Paragraph("para"));
	
	% S-S distances
	add(rpt, Heading1("Cys-Cys távolságok"));
	res = unique([this.DistancesS.res1; this.DistancesS.res2]);
	chain = unique(this.DistancesS.chain);
	for iRes1 = 1:length(res)
		add(rpt, Heading2(strcat(num2str(res(iRes1)), "C")));
		% Group by chains
		for iChain = 1:length(chain)
			add(rpt, Heading3(strcat(num2str(res(iRes1)), "C-xC (chain ", chain(iChain), ")")));
			figure('units', 'centimeters', 'position', [0 0 15 9]);
			leg = string.empty;
			for iRes2 = 1:length(res)
				if iRes1 ~= iRes2
					ind = find((this.DistancesS.res1 == res(iRes1) | this.DistancesS.res2 == res(iRes1)) & (this.DistancesS.res1 == res(iRes2) | this.DistancesS.res2 == res(iRes2)) & this.DistancesS.chain == chain(iChain), 1);
					switch res(iRes2)
						case 54
							color = [1 0 0];
						case 61
							color = [0.5 0 0];
						case 65
							color = [1 0.5 0];
						case 187
							color = [0 1 0];
						case 192
							color = [0 0.5 0];
						case 198
							color = [0.5 1 0];
					end
					plot(this.Time, cell2mat(this.DistancesS.dist(ind)), '.', 'Color', color);
					hold on;
					xlabel('Time (ns)');
					ylabel('S-S distance (A)');
					ylim([1  6]);
					leg = [leg, strcat(num2str(res(iRes1)), 'C-', num2str(res(iRes2)), 'C')];
				end
			end
			legend(leg, 'Location', 'northoutside', 'NumColumns', 5);
			fig = Figure;
			fig.SnapshotFormat = "png";
			add(rpt, fig);
			close(gcf);
		end
		
		% Group by Cys pairs
		for iRes2 = 1:length(res)
			if iRes1 ~= iRes2
				add(rpt, Heading3(strcat(num2str(res(iRes1)), "C-", num2str(res(iRes2)), "C (all chains)")));
				figure('units', 'centimeters', 'position', [0 0 15 9]);
				leg = string.empty;
				for iChain = 1:length(chain)
					ind = find((this.DistancesS.res1 == res(iRes1) | this.DistancesS.res2 == res(iRes1)) & (this.DistancesS.res1 == res(iRes2) | this.DistancesS.res2 == res(iRes2)) & this.DistancesS.chain == chain(iChain), 1);
					plot(this.Time, cell2mat(this.DistancesS.dist(ind)), '.');
					hold on;
					xlabel('Time (ns)');
					ylabel('S-S distance (A)');
					ylim([1 6]);
					leg = [leg, strcat("Chain ", chain(iChain))];
				end
				legend(leg, 'Location', 'northoutside', 'NumColumns', 6);
				fig = Figure;
				fig.SnapshotFormat = "png";
				add(rpt, fig);
				close(gcf);
			end
		end
	end
	
	% Calpha distances
	add(rpt, Heading1("EL1-EL2 távolságok"));
	chain = unique([this.DistancesCA.chain1; this.DistancesCA.chain1]);
	add(rpt, Heading2("58Q-193P (all chains)"));
	figure('units', 'centimeters', 'position', [0 0 15 9]);
	leg = string.empty;
	for iChain = 1:length(chain)
		ind = find((this.DistancesCA.res1 == 58 | this.DistancesCA.res2 == 58) & (this.DistancesCA.res1 == 193 | this.DistancesCA.res2 == 193) & (this.DistancesCA.chain1 == chain(iChain) | this.DistancesCA.chain2 == chain(iChain)), 1);
		plot(this.Time, cell2mat(this.DistancesCA.dist(ind)), '.');
		hold on;
		xlabel('Time (ns)');
		ylabel('C_\alpha-C_\alpha distance (A)');
		leg = [leg, strcat("Chain ", chain(iChain))];
	end
	legend(leg, 'Location', 'northoutside', 'NumColumns', 6);
	fig = Figure;
	fig.SnapshotFormat = "png";
	add(rpt, fig);
	close(gcf);
		
	if max(this.ResiduePositionCA.chain == "G") == 1
		add(rpt, Heading1("Transz-junkcionális távolságok"));
		chains = {'A', 'J'; 'B', 'I'; 'C', 'H'; 'D', 'G'; 'E', 'L'; 'F', 'K'};
		for iChain = 1:size(chains, 1)
			add(rpt, Heading2(strcat("Chains ", chains(iChain, 1), "-", chains(iChain, 2))));
			figure('units', 'centimeters', 'position', [0 0 15 9]);

			posCA = this.ResiduePositionCA(this.ResiduePositionCA.chain == chains(iChain, 1),:);
			z = reshape(posCA.z, max(posCA.frame) + 1, length(unique(posCA.residue)), []);
			plot((1:size(z, 1)) ./ this.FrameInterval, max(z, [], 2));
			hold on;

			posCA = this.ResiduePositionCA(this.ResiduePositionCA.chain == chains(iChain, 2),:);
			z = reshape(posCA.z, max(posCA.frame) + 1, length(unique(posCA.residue)), []);
			plot((1:size(z, 1)) ./ this.FrameInterval, min(z, [], 2));

			xlabel('Time (ns)');
			ylabel('z coordinate');
			leg = {strcat("max (", chains(iChain, 1), " chain)"), strcat("min (", chains(iChain, 2), " chain)")};
			legend(leg, 'Location', 'northoutside', 'NumColumns', 3);

			fig = Figure;
			fig.SnapshotFormat = "png";
			add(rpt, fig);
			close(gcf);
		end
	end
		
	close(rpt);
	rptview(rpt);
end