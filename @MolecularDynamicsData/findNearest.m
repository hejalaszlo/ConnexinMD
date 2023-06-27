function result = findNearest(this, resid)
	result = table.empty;
	chain = unique(this.ResiduePositionCA.chain);
	res = unique(this.ResiduePositionCA.residue);
	frame = unique(this.ResiduePositionCA.frame);
	for iChain1 = 1:length(chain)
		for iFrame = min(frame):max(frame)
			x1 = this.ResiduePositionCA.x(this.ResiduePositionCA.chain == chain(iChain1) & this.ResiduePositionCA.frame == iFrame & this.ResiduePositionCA.residue == resid);
			y1 = this.ResiduePositionCA.y(this.ResiduePositionCA.chain == chain(iChain1) & this.ResiduePositionCA.frame == iFrame & this.ResiduePositionCA.residue == resid);
			z1 = this.ResiduePositionCA.z(this.ResiduePositionCA.chain == chain(iChain1) & this.ResiduePositionCA.frame == iFrame & this.ResiduePositionCA.residue == resid);
			
			x2 = this.ResiduePositionCA.x(this.ResiduePositionCA.chain == chain(iChain1) & this.ResiduePositionCA.frame == iFrame & abs(double(this.ResiduePositionCA.residue) - resid) > 50);
			y2 = this.ResiduePositionCA.y(this.ResiduePositionCA.chain == chain(iChain1) & this.ResiduePositionCA.frame == iFrame & abs(double(this.ResiduePositionCA.residue) - resid) > 50);
			z2 = this.ResiduePositionCA.z(this.ResiduePositionCA.chain == chain(iChain1) & this.ResiduePositionCA.frame == iFrame & abs(double(this.ResiduePositionCA.residue) - resid) > 50);
			res2 = this.ResiduePositionCA.residue(this.ResiduePositionCA.chain == chain(iChain1) & this.ResiduePositionCA.frame == iFrame & abs(double(this.ResiduePositionCA.residue) - resid) > 50);
			
			dist = sqrt((x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2);
			minid = find(dist == min(dist), 1);
			result = [result; table(chain(iChain1), iFrame, dist(minid), res2(minid))];
		end
	end
	result.Properties.VariableNames = {'chain' 'frame' 'dist' 'res'};
end