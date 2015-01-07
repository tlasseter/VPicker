package com.Sts.Actions.Wizards.VolumeStimulated;

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

public class StsSelectSubVolumesPanel extends StsJPanel
{
    private StsVolumeStimulatedWizard wizard;
    private StsSelectSubVolumes wizardStep;

    private StsSubVolume[] availableSubVolumes = null;    
    private StsModel model = null;

    StsJPanel beanPanel = new StsJPanel();
    StsListFieldBean subVolumeListBean;
    
    public StsSelectSubVolumesPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsVolumeStimulatedWizard)wizard;
    	this.wizardStep = (StsSelectSubVolumes)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
        if(subVolumeClass != null && subVolumeClass.getIsVisible())
        {
            availableSubVolumes = subVolumeClass.getSubVolumes();
    	    subVolumeListBean = new StsListFieldBean(wizard, "subVolumeObjects", "SubVolumes: ", availableSubVolumes);
    	    subVolumeListBean.setSelectedAll();
        }

        gbc.fill = gbc.BOTH;
    	add(subVolumeListBean);    	
    }
    
    public Object[] getSelectedObjects()
    {
    	return subVolumeListBean.getSelectedObjects();
    }
}
