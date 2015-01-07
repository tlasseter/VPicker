package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineSubVolumePanel extends StsJPanel
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;
    
    protected StsGroupBox boxPanel = new StsGroupBox("Define");
//	protected StsComboBoxFieldBean domainComboBean = new StsComboBoxFieldBean();
	protected StsComboBoxFieldBean typesComboBean = new StsComboBoxFieldBean();
	protected StsStringFieldBean cvNameBean = new StsStringFieldBean();
	
	protected byte subVolumeType = 0;
	protected String svName = "None";
//	protected byte domainType = 0;
	
    public StsDefineSubVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
        /*
        domainType = model.getProject().getZDomain();
        domainComboBean.setSelectedIndex(domainType);        
        if(model.getProject().getZDomainSupported() != StsParameters.TD_TIME_DEPTH)
        	domainComboBean.setEditable(false);
        */
    }

    void jbInit() throws Exception
    {
    	typesComboBean.initialize(this, "subVolumeTypeString", "Type", StsSubVolume.svTypes);
    	typesComboBean.setToolTipText("Select the type of subvolume you wish to create.");

        // domainComboBean.initialize(this, "DomainString", "Domain:", StsParameters.TD_STRINGS);
    	// domainComboBean.setToolTipText("Select the domain which you will be defining the subvolume for.");
    	cvNameBean.initialize(this, "subVolumeName", true, "Name:");
    	
    	boxPanel.gbc.fill = gbc.HORIZONTAL;
    	boxPanel.addEndRow(cvNameBean);
    	boxPanel.addEndRow(typesComboBean);
    	// boxPanel.addEndRow(domainComboBean);
    	
    	gbc.fill = gbc.HORIZONTAL;
    	gbc.anchor = gbc.NORTH;
    	add(boxPanel);
    }

    public byte getSubVolumeType()
    {
        return subVolumeType;
    }
    public String getSubVolumeTypeString()
    {
        return StsSubVolume.svTypes[subVolumeType];
    }    
    public void setSubVolumeTypeString(String typeString)
    {
    	subVolumeType = StsSubVolume.getSubVolumeTypeFromString(typeString);
    }
    /*
    public byte getDomain()
    {
        return domainType;
    }
    public String getDomainString()
    {
        return StsParameters.TD_STRINGS[domainType];
    }    
    public void setDomainString(String domainString)
    {
    	domainType = StsParameters.getZDomainFromString(domainString);
    }
    */
    public String getSubVolumeName()
    {
        return svName;
    }
    public void setSubVolumeName(String subName)
    {
        svName = subName;
    }    
}
