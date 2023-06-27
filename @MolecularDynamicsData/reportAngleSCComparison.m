function reportAngleSCComparison(this)
	import mlreportgen.report.* 
	import mlreportgen.dom.* 
	
	rpt = Report(strcat("Torsion angles vs. SCs, ", inputname(1)), "docx");

	add(rpt, Heading1("Összefoglalás"));
	add(rpt, Paragraph("para"));
	
	% Torsion angles, SC dynamics
	add(rpt, Heading1("Prolin torziós szögek, SC dinamika"));
	res = unique(this.TorsionAngles.residue);
	chain = unique(this.TorsionAngles.chain);
	% Group by chain
% 	for iChain = 1:length(chain)
% 		add(rpt, Heading2(strcat("Chain ", chain(iChain))));
% 		figure('units', 'centimeters', 'position', [0 0 15 15]);
% 		
% 		% SCs
% 		subplot(3,1,1);
% 		stabCenters = this.StabCentersExtracellFrequent(strcmp(this.StabCentersExtracellFrequent.Chain1, chain(iChain)) == 1 | strcmp(this.StabCentersExtracellFrequent.Chain2, chain(iChain)) == 1, :);
% 		resovertime = false(height(stabCenters), this.FrameNum);
% 		for iSC = 1:height(stabCenters)
% 			resovertime(iSC, stabCenters.Frames{iSC}) = true;
% 		end
% 		interchain = strcmp(stabCenters.Chain1, stabCenters.Chain2) == 0;
% 		interchainNum = sum(interchain);
% 
% 		if interchainNum > 0
% 			imagesc(this.Time, 1:height(stabCenters(interchain,:)), resovertime(interchain,:));
% 			title('Stab. centers, interchain, occurrence > 1 %');
% 			cmapb = ones(2,3);
% 			cmapb(2,:) = [0 0 1];
% 			colormap(cmapb);
% 			xlabel('Time (ns)');
% 			set(gca, 'Ytick', 1:height(stabCenters(interchain,:)), 'YTickLabel', strcat(stabCenters.Chain1(interchain), num2str(stabCenters.Pos1(interchain)), stabCenters.Res1(interchain), '-', stabCenters.Chain2(interchain), num2str(stabCenters.Pos2(interchain)), stabCenters.Res2(interchain)));
% 		end
% 
% 		% Torsion angles
% 		subplot(3,1,[2 3]);
% 		for iRes = 1:length(res)
% 			ind = this.TorsionAngles.chain == chain(iChain) & this.TorsionAngles.residue == res(iRes);
% 			plot(double(this.TorsionAngles.frame(ind))/10, this.TorsionAngles.psi(ind), '.');
% 			hold on;
% 		end
% 		xlabel('Time (ns)');
% 		ylabel('Torsion angle (degree)');
% 		leg = {'Pro59 psi', 'Pro191 psi', 'Pro193 psi'};
% 		legend(leg, 'Location', 'northoutside', 'NumColumns', 3);
% 		fig = Figure;
% 		fig.SnapshotFormat = "png";
% 		add(rpt, fig);
% 		close(gcf);
% 	end

	% Group by chain pairs
	chains = {'A', 'J'; 'B', 'I'; 'C', 'H'; 'D', 'G'; 'E', 'L'; 'F', 'K'};
	for iChain = 1:size(chains, 1)
		add(rpt, Heading2(strcat("Chains ", chains(iChain, 1), "-", chains(iChain, 2))));
		figure('units', 'centimeters', 'position', [0 0 15 18]);
		
		% SCs
		subplot(5,1,1);
		stabCenters = this.StabCentersTransjunctional((strcmp(this.StabCentersTransjunctional.Chain1, chains(iChain, 1)) == 1 & strcmp(this.StabCentersTransjunctional.Chain2, chains(iChain, 2)) == 1) | (strcmp(this.StabCentersTransjunctional.Chain1, chains(iChain, 2)) == 1 & strcmp(this.StabCentersTransjunctional.Chain2, chains(iChain, 1)) == 1), :);
		resovertime = false(height(stabCenters), this.FrameNum);
		for iSC = 1:height(stabCenters)
			resovertime(iSC, stabCenters.Frames{iSC}) = true;
		end
		imagesc(this.Time, 1:height(stabCenters), resovertime);
		title('Stab. centers, transjunctional');
		cmapb = ones(2,3);
		cmapb(2,:) = [0 0 1];
		colormap(cmapb);
		xlabel('Time (ns)');
		set(gca, 'Ytick', 1:height(stabCenters), 'YTickLabel', strcat(stabCenters.Chain1, num2str(stabCenters.Pos1), stabCenters.Res1, '-', stabCenters.Chain2, num2str(stabCenters.Pos2), stabCenters.Res2));

		% Torsion angles, chain 1
		subplot(5,1,[2 3]);
		for iRes = 1:length(res)
			ind = this.TorsionAngles.chain == chains(iChain, 1) & this.TorsionAngles.residue == res(iRes);
			plot(double(this.TorsionAngles.frame(ind))/10, this.TorsionAngles.psi(ind), '.');
			hold on;
		end
		ylabel('Psi torsion angle (degree)');
		leg = {strcat("Chain ", chains(iChain, 1), " 59P"), strcat("Chain ", chains(iChain, 1), " 191P"), strcat("Chain ", chains(iChain, 1), " 193P")};
		legend(leg, 'Location', 'northoutside', 'NumColumns', 3);
		
		% Torsion angles, chain 2
		subplot(5,1,[4 5]);
		for iRes = 1:length(res)
			ind = this.TorsionAngles.chain == chains(iChain, 2) & this.TorsionAngles.residue == res(iRes);
			plot(double(this.TorsionAngles.frame(ind))/10, this.TorsionAngles.psi(ind), '.');
			hold on;
		end
		xlabel('Time (ns)');
		ylabel('Psi torsion angle (degree)');
		leg = {strcat("Chain ", chains(iChain, 2), " 59P"), strcat("Chain ", chains(iChain, 2), " 191P"), strcat("Chain ", chains(iChain, 2), " 193P")};
		legend(leg, 'Location', 'northoutside', 'NumColumns', 3);

		fig = Figure;
		fig.SnapshotFormat = "png";
		add(rpt, fig);
		close(gcf);
	end

	close(rpt);
	rptview(rpt);
end