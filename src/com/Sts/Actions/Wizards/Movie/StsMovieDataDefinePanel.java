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

public class StsMovieDataDefinePanel extends JPanel
{
    private StsMovieWizard wizard;
    private StsMovieDataDefine wizardStep;

    JRadioButton inlineRadio = new JRadioButton();
    JPanel jPanel1 = new JPanel();
    JRadioButton xLineRadio = new JRadioButton();
    JRadioButton zRadio = new JRadioButton();
    ButtonGroup orientGrp = new ButtonGroup();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
  GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsMovieDataDefinePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMovieWizard)wizard;
        this.wizardStep = (StsMovieDataDefine)wizardStep;

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

        inlineRadio.setSelected(true);
        jPanel1.add(inlineRadio,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 5, 0, 5), 0, 0));
        jPanel1.add(xLineRadio,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 5, 0, 5), 0, 0));
        jPanel1.add(zRadio,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 5, 8, 5), 0, 0));


        this.add(jPanel1,   new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 10, 7, 14), 64, -3));

        inlineRadio.setText("In Line");

        this.setMinimumSize(new Dimension(300, 300));
        jPanel1.setBorder(BorderFactory.createEtchedBorder());

        xLineRadio.setText("Cross Line");
        zRadio.setText("Time / Depth");
        orientGrp.add(inlineRadio);
        orientGrp.add(xLineRadio);
        orientGrp.add(zRadio);
    }

    public int getDirection()
    {
        if(inlineRadio.isSelected()) return StsMovie.INLINE;
        if(xLineRadio.isSelected()) return StsMovie.XLINE;
        if(zRadio.isSelected()) return StsMovie.ZDIR;
        return (int)-1;
    }
}
