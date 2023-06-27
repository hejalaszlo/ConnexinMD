function reportAngles(this)
	import mlreportgen.report.* 
	import mlreportgen.dom.* 
	
	rpt = Report(strcat("Torsion angles, ", inputname(1)), "docx");

	add(rpt, Heading1("Összefoglalás"));
	add(rpt, Paragraph("para"));
	
	% S-S distances
	add(rpt, Heading1("Prolin torziós szögek"));
	res = unique(this.TorsionAngles.residue);
	chain = unique(this.TorsionAngles.chain);
	for iRes = 1:length(res)
		add(rpt, Heading2(strcat("Pro", num2str(res(iRes)))));
		% Group by chain
		for iChain = 1:length(chain)
			add(rpt, Heading3(strcat("Chain ", chain(iChain))));
			figure('units', 'centimeters', 'position', [0 0 15 9]);
			ind = this.TorsionAngles.chain == chain(iChain) & this.TorsionAngles.residue == res(iRes);
			plot(double(this.TorsionAngles.frame(ind))/10, this.TorsionAngles.phi(ind), '.');
			hold on;
			plot(double(this.TorsionAngles.frame(ind))/10, this.TorsionAngles.psi(ind), '.');
			plot(double(this.TorsionAngles.frame(ind))/10, this.TorsionAngles.omegaprev(ind), '.');
			xlabel('Time (ns)');
			ylabel('Torsion angle (degree)');
			ylim([-200 200]);
			
			leg = {'Phi', 'Psi', 'Omega (prev)'};
			legend(leg, 'Location', 'northoutside', 'NumColumns', 3);
			fig = Figure;
			fig.SnapshotFormat = "png";
			add(rpt, fig);
			close(gcf);
		end
			
		% Group by angle
		add(rpt, Heading3("Phi"));
		figure('units', 'centimeters', 'position', [0 0 15 9]);
		ind = this.TorsionAngles.residue == res(iRes);
		gscatter(double(this.TorsionAngles.frame(ind))/10, this.TorsionAngles.phi(ind), this.TorsionAngles.chain(ind));
		xlabel('Time (ns)');
		ylabel('Torsion angle (degree)');
		ylim([-200 200]);
		legend('Location', 'northoutside', 'NumColumns', length(chain));
		fig = Figure;
		fig.SnapshotFormat = "png";
		add(rpt, fig);
		close(gcf);
		
		add(rpt, Heading3("Psi"));
		figure('units', 'centimeters', 'position', [0 0 15 9]);
		ind = this.TorsionAngles.residue == res(iRes);
		gscatter(double(this.TorsionAngles.frame(ind))/10, this.TorsionAngles.psi(ind), this.TorsionAngles.chain(ind));
		xlabel('Time (ns)');
		ylabel('Torsion angle (degree)');
		ylim([-200 200]);
		legend('Location', 'northoutside', 'NumColumns', length(chain));
		fig = Figure;
		fig.SnapshotFormat = "png";
		add(rpt, fig);
		close(gcf);
		
		add(rpt, Heading3("Omega (prev)"));
		figure('units', 'centimeters', 'position', [0 0 15 9]);
		ind = this.TorsionAngles.residue == res(iRes);
		gscatter(double(this.TorsionAngles.frame(ind))/10, this.TorsionAngles.omegaprev(ind), this.TorsionAngles.chain(ind));
		xlabel('Time (ns)');
		ylabel('Torsion angle (degree)');
		ylim([-200 200]);
		legend('Location', 'northoutside', 'NumColumns', length(chain));
		fig = Figure;
		fig.SnapshotFormat = "png";
		add(rpt, fig);
		close(gcf);
	end
		
	close(rpt);
	rptview(rpt);
end