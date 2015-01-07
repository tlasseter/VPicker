package com.Sts.Actions.Wizards.Movie;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;

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

public class StsMovieVolumeDefinePanel extends JPanel
{
    private StsMovieWizard wizard;
    private StsMovieVolumeDefine wizardStep;

    JPanel jPanel1 = new JPanel();
    JRadioButton pointSetRadio = new JRadioButton();
    JRadioButton voxelRadio = new JRadioButton();
    ButtonGroup volumeGrp = new ButtonGroup();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsMovieVolumeDefinePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMovieWizard)wizard;
        this.wizardStep = (StsMovieVolumeDefine)wizardStep;

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

    }

    void jbInit() throws Exception
    {
         this.setLayout(gridBagLayout1);
        jPanel1.setLayout(gridBagLayout2);

        pointSetRadio.setSelected(true);
        jPanel1.add(pointSetRadio,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 5, 0, 5), 0, 0));
        jPanel1.add(voxelRadio,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 5, 8, 5), 0, 0));


        this.add(jPanel1,   new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 10, 7, 14), 64, -3));

        this.setMinimumSize(new Dimension(300, 300));
        jPanel1.setBorder(BorderFactory.createEtchedBorder());

        pointSetRadio.setText("Point Set");
        voxelRadio.setText("Voxel");
        volumeGrp.add(voxelRadio);
        volumeGrp.add(pointSetRadio);
    }

    public byte getVolumeType()
    {
        if(voxelRadio.isSelected()) return StsMovie.VOXEL;
        if(pointSetRadio.isSelected()) return StsMovie.POINTSET;
        return (int)-1;
    }
}
