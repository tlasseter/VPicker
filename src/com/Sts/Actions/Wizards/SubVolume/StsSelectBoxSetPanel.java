package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectBoxSetPanel extends JPanel
{
    private StsWizard wizard;
    private StsSelectBoxSet selectBoxSet;

    StsListFieldBean boxSetListBean = new StsListFieldBean();
    StsButton newBoxSetButton = new StsButton();
    JLabel listLabel = new JLabel("Available Box Set SubVolumes");

    public StsSelectBoxSetPanel(StsWizard wizard, StsSelectBoxSet selectBoxStep)
    {
        this.wizard = wizard;
        this.selectBoxSet = selectBoxStep;
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
        StsModel model = wizard.getModel();
        StsBoxSetSubVolumeClass boxSetClass = (StsBoxSetSubVolumeClass)model.getStsClass(StsBoxSetSubVolume.class);
        StsBoxSetSubVolume[] boxSets = (StsBoxSetSubVolume[])boxSetClass.getCastObjectList(StsBoxSetSubVolume.class);
		boxSetListBean.initialize(selectBoxSet, "boxSet", false);
		boxSetListBean.setListItems(boxSets);
		StsBoxSetSubVolume currentBoxSet = (StsBoxSetSubVolume)boxSetClass.getCurrentObject();
       if(currentBoxSet != null) boxSetListBean.setSelectedValue(currentBoxSet);
    }

    void jbInit() throws Exception
    {
        setLayout(new GridBagLayout());

        boxSetListBean.initialize(selectBoxSet, "boxSet", "", null);
        boxSetListBean.setBorder(BorderFactory.createEtchedBorder());

        this.add(listLabel, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 5, 2, 0), 0, 0));
        this.add(boxSetListBean,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(0, 0, 5, 0), 150, 0));

        newBoxSetButton.initialize("New Box Set SubVolume", "Create a set of hexahedral subvolumes.", selectBoxSet, "createBoxSet", null);

        this.add(newBoxSetButton,   new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 100, 0));
    }
}
