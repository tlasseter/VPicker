package com.Sts.Actions.Wizards.Movie;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsMovieTypePanel extends StsJPanel
{
    private StsMovieWizard wizard;
    private StsMovieType wizardStep;

    StsGroupBox nameGroup = new StsGroupBox();
    StsStringFieldBean movieNameBean = new StsStringFieldBean();
    String movieName = "MovieName";

    StsGroupBox orientBox = new StsGroupBox("Orientation");
    JRadioButton inlineBean;
    JRadioButton xlineBean;
    JRadioButton sliceBean;
    ButtonGroup orientGrp = new ButtonGroup();

    public StsMovieTypePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMovieWizard)wizard;
        this.wizardStep = (StsMovieType)wizardStep;

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
        gbc.fill = gbc.HORIZONTAL;
        movieNameBean = new StsStringFieldBean(this, "MovieName", true, "Name:");
        nameGroup.gbc.fill = gbc.HORIZONTAL;
        nameGroup.add(movieNameBean);
        inlineBean = new JRadioButton("Inline", false);
        xlineBean = new JRadioButton("Crossline", false);
        sliceBean = new JRadioButton("Slice", true);
        orientGrp.add(inlineBean);
        orientGrp.add(xlineBean);
        orientGrp.add(sliceBean);
        orientBox.gbc.fill = gbc.HORIZONTAL;
        orientBox.addEndRow(inlineBean);
        orientBox.addEndRow(xlineBean);
        orientBox.addEndRow(sliceBean);
        addEndRow(nameGroup);
        addEndRow(orientBox);
    }

    public String getMovieName()
    {
        return movieName;
    }

    public void setMovieName(String name)
    {
        movieName = name;
    }

    public int getDirection()
    {
        if(inlineBean.isSelected()) return StsMovie.INLINE;
        if(xlineBean.isSelected()) return StsMovie.XLINE;
        return StsMovie.ZDIR;
    }
}
