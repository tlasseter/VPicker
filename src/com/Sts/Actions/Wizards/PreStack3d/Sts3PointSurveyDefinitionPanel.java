package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 *
 * @author unascribed
 * @version 1.0
 */

public class Sts3PointSurveyDefinitionPanel extends StsJPanel
{
    StsSeismicBoundingBox volume;
    StsGroupBox knownBox = new StsGroupBox("Known Points");
    StsGroupBox graphicBox = new StsGroupBox("Survey Graphic");

    StsSegyVolume[] vols = null;

    static final float nullValue = StsParameters.nullValue;
    static final double doubleNullValue = StsParameters.doubleNullValue;

    public Sts3PointSurveyDefinitionPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        try
        {
            buildPanel();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize(StsSeismicBoundingBox volume)
    {
        this.volume = volume;
        volume.initializeKnownPoints();
        knownBox.setPanelObject(volume);
    }

    private void buildPanel()
    {
        StsDoubleFieldBean point1X = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "point1X", "Point 1: X");
        StsDoubleFieldBean point1Y = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "point1Y", "Y");
        StsFloatFieldBean point1Inline = new StsFloatFieldBean(StsSeismicBoundingBox.class, "point1Inline", "Inline");
        StsFloatFieldBean point1Xline = new StsFloatFieldBean(StsSeismicBoundingBox.class, "point1Xline", "Xline");

        StsDoubleFieldBean point2X = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "point2X", "Point 2: X");
        StsDoubleFieldBean point2Y = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "point2Y", "Y");
        StsFloatFieldBean point2Inline = new StsFloatFieldBean(StsSeismicBoundingBox.class, "point2Inline", "Inline");
        StsFloatFieldBean point2Xline = new StsFloatFieldBean(StsSeismicBoundingBox.class, "point2Xline", "Xline");

        StsDoubleFieldBean point3X = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "point3X", "Point 3: X");
        StsDoubleFieldBean point3Y = new StsDoubleFieldBean(StsSeismicBoundingBox.class, "point3Y", "Y");
        StsFloatFieldBean point3Inline = new StsFloatFieldBean(StsSeismicBoundingBox.class, "point3Inline", "Inline");
        StsFloatFieldBean point3Xline = new StsFloatFieldBean(StsSeismicBoundingBox.class, "point3Xline", "Xline");

        knownBox.addToRow(point1X);
        knownBox.addToRow(point1Y);
        knownBox.addToRow(point1Inline);
        knownBox.addEndRow(point1Xline);
        knownBox.addToRow(point2X);
        knownBox.addToRow(point2Y);
        knownBox.addToRow(point2Inline);
        knownBox.addEndRow(point2Xline);
        knownBox.addToRow(point3X);
        knownBox.addToRow(point3Y);
        knownBox.addToRow(point3Inline);
        knownBox.addEndRow(point3Xline);
        add(knownBox);
 /*
        graphicBox.setBackground(Color.white);
        graphicBox.setFont(new Font("Dialog", 1, 11));
        graphicBox.setPreferredSize(new Dimension(500, 250));
        add(graphicBox);
 */
    }
}
