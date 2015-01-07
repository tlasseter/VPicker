package com.Sts.Actions.Wizards.PreStackExport;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.*;

public class StsExportHandvel2d extends StsExportHandvel 
{

	public StsExportHandvel2d(StsWizard wizard) 
	{
		super(wizard);
	}

	@Override
	protected StsPreStackHandVelExportPanel createHandVelPanel() 
	{
		return new StsPreStackHandVelExportPanel2d(model, velocityModel);
	}

}
