function reportSC(this)
	import mlreportgen.report.* 
	import mlreportgen.dom.* 
	
	rpt = Report(strcat("SC, ", inputname(1)), "docx");

	add(rpt, Heading1("Összefoglalás"));
	add(rpt, Paragraph("Olyan SC-k, amik a futásidő legalább 2 %-ában jelen vannak."));
	
	% SC-k száma
	add(rpt, Heading1("SC-k száma"));

	stabCenters = this.StabCenters;
	resovertime = false(height(stabCenters), this.FrameNum);
	for iRes = 1:height(stabCenters)
		resovertime(iRes, stabCenters.Frames{iRes}) = true;
	end
	
	figure;
	plot(this.Time, sum(resovertime)');
	hold on;

	stabCenters = this.StabCentersFrequent;
	resovertime = false(height(stabCenters), this.FrameNum);
	for iRes = 1:height(stabCenters)
		resovertime(iRes, stabCenters.Frames{iRes}) = true;
	end
	
	plot(this.Time, sum(resovertime)');
	xlabel('Time (ns)');
	legend({'Minden SC', 'Stabil SC-k (> 2 %)'});

	fig = Figure;
	fig.SnapshotFormat = "png";
	add(rpt, fig);
	close(gcf);

	% Elhelyezkedés szerint
	add(rpt, Heading1("Elhelyezkedés szerint"));

	% Trans-junctional SCs
	add(rpt, Heading2("Trans-junctional"));
	ind = (ismember(stabCenters.Chain1, {'A', 'B', 'C', 'D', 'E', 'F'}) & ismember(stabCenters.Chain2, {'G', 'H', 'I', 'J', 'K', 'L'})) | (ismember(stabCenters.Chain2, {'A', 'B', 'C', 'D', 'E', 'F'}) & ismember(stabCenters.Chain1, {'G', 'H', 'I', 'J', 'K', 'L'}));
	if sum(ind) > 0
		figure;
		imagesc(this.Time, 1:height(stabCenters(ind,:)), resovertime(ind,:));
		title('Stab. centers, transjunction');
		cmapb = ones(2,3);
		cmapb(2,:) = [0 0 1];
		colormap(cmapb);
		xlabel('Time (ns)');
		set(gca, 'Ytick', 1:height(stabCenters(ind,:)), 'YTickLabel', strcat(stabCenters.Chain1(ind), num2str(stabCenters.Pos1(ind)), stabCenters.Res1(ind), '-', stabCenters.Chain2(ind), num2str(stabCenters.Pos2(ind)), stabCenters.Res2(ind)));
		
		fig = Figure;
		fig.SnapshotFormat = "png";
		add(rpt, fig);
		close(gcf);
	end

	% Inter-subunit SCs
	add(rpt, Heading2("Inter-subunit"));
	ind = ~ind & strcmp(stabCenters.Chain1, stabCenters.Chain2) == 0;
	if sum(ind) > 0	
		figure;
		imagesc(this.Time, 1:height(stabCenters(ind,:)), resovertime(ind,:));
		title('Stab. centers, intersubunit');
		cmapb = ones(2,3);
		cmapb(2,:) = [0 1 0];
		colormap(cmapb);
		xlabel('Time (ns)');
		set(gca, 'Ytick', 1:height(stabCenters(ind,:)), 'YTickLabel', strcat(stabCenters.Chain1(ind), num2str(stabCenters.Pos1(ind)), stabCenters.Res1(ind), '-', stabCenters.Chain2(ind), num2str(stabCenters.Pos2(ind)), stabCenters.Res2(ind)));
			
		fig = Figure;
		fig.SnapshotFormat = "png";
		add(rpt, fig);
		close(gcf);
	end

	% Inter-loop SCs
	add(rpt, Heading2("Inter-loop"));
	ind = strcmp(stabCenters.Chain1, stabCenters.Chain2) == 1 & abs(stabCenters.Pos1 - stabCenters.Pos2) > 50;
	if sum(ind) > 0	
		sc = stabCenters(ind,:);
		rov = resovertime(ind,:);
		for i = 1:30:height(sc)
			figure('units', 'centimeters', 'position', [0 0 15 9]);
			sci = sc(i:min(i+29,height(sc)),:);
			rovi = rov(i:min(i+29,height(sc)),:);
			imagesc(this.Time, 1:height(sci), rovi);
			title('Stab. centers, interloop');
			cmapr = ones(2,3);
			cmapr(2,:) = [1 0 0];
			colormap(cmapr);
			xlabel('Time (ns)');
			set(gca, 'Ytick', 1:height(sci), 'YTickLabel', strcat(sci.Chain1, num2str(sci.Pos1), sci.Res1, '-', sci.Chain2, num2str(sci.Pos2), sci.Res2));
			
			fig = Figure;
			fig.SnapshotFormat = "png";
			add(rpt, fig);
			close(gcf);
		end
	end
		
	% CYD motívum
	add(rpt, Heading2("CYD motívum"));
	add(rpt, Paragraph("A 65C-66Y-67D aminosavak egyikét tartalmazó SC-k."));
	ind = stabCenters.Pos1 == 65 | stabCenters.Pos1 == 66 | stabCenters.Pos1 == 67 | stabCenters.Pos2 == 65 | stabCenters.Pos2 == 62 | stabCenters.Pos2 == 67;
	if sum(ind) > 0	
		figure;
		imagesc(this.Time, 1:height(stabCenters(ind,:)), resovertime(ind,:));
		title('Stab. centers, CYD motívum');
		cmapb = ones(2,3);
		cmapb(2,:) = [0 0 1];
		colormap(cmapb);
		xlabel('Time (ns)');
		set(gca, 'Ytick', 1:height(stabCenters(ind,:)), 'YTickLabel', strcat(stabCenters.Chain1(ind), num2str(stabCenters.Pos1(ind)), stabCenters.Res1(ind), '-', stabCenters.Chain2(ind), num2str(stabCenters.Pos2(ind)), stabCenters.Res2(ind)));
			
		fig = Figure;
		fig.SnapshotFormat = "png";
		add(rpt, fig);
		close(gcf);
	end
		
	% Aminosav típusa szerint
	add(rpt, Heading1("Aminosav típusa szerint"));

	% Lysine-cysteine SCs
	add(rpt, Heading2("Lysine-Cysteine"));
	ind = (stabCenters.Res1 == "K" & stabCenters.Res2 == "C") | (stabCenters.Res1 == "C" & stabCenters.Res2 == "K");
	if sum(ind) > 0
		figure;
		imagesc(this.Time, 1:height(stabCenters(ind,:)), resovertime(ind,:));
		title('Stab. centers, Lysine-Cysteine');
		cmapb = ones(2,3);
		cmapb(2,:) = [0 0 1];
		colormap(cmapb);
		xlabel('Time (ns)');
		set(gca, 'Ytick', 1:height(stabCenters(ind,:)), 'YTickLabel', strcat(stabCenters.Chain1(ind), num2str(stabCenters.Pos1(ind)), stabCenters.Res1(ind), '-', stabCenters.Chain2(ind), num2str(stabCenters.Pos2(ind)), stabCenters.Res2(ind)));
		
		fig = Figure;
		fig.SnapshotFormat = "png";
		add(rpt, fig);
		close(gcf);
	end
	
	close(rpt);
	rptview(rpt);
end