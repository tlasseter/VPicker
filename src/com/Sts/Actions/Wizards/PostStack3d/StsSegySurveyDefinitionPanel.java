package com.Sts.Actions.Wizards.PostStack3d;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSegySurveyDefinitionPanel extends StsJPanel
{
        StsGroupBox knownBox = new StsGroupBox("Known Points");
        StsGroupBox graphicBox = new StsGroupBox("Survey Graphic");
        StsGroupBox buttonBox = new StsGroupBox();

        private StsSeismicWizard wizard;
        private StsWizardStep wizardStep;
        private JDialog parent;

        StsSeismicBoundingBox[] vols = null;

        static final float nullValue = StsParameters.nullValue;
        static final double doubleNullValue = StsParameters.doubleNullValue;

        public StsSegySurveyDefinitionPanel(StsSeismicWizard wizard, StsWizardStep wizardStep, JDialog dialog)
        {
            this.wizard = wizard;
            this.wizardStep = wizardStep;
            this.parent = dialog;
            try
            {
                vols = wizard.getSegyVolumes();
                initialize();
                if(vols != null && vols[0] != null)
                {
                    vols[0].initializeKnownPoints();
                    knownBox.setPanelObject(vols[0]);
                }
            }
            catch(Exception ex) {
                ex.printStackTrace();
            }
        }

        private void initialize()
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

            StsButton acceptBtn = new StsButton("Accept", "Recompute geometry based on this trace info.", this, "accept");
            StsButton cancelBtn = new StsButton("Cancel", "Ignore this trace info.", this, "cancel");

            knownBox.addToRow(point1X); knownBox.addToRow(point1Y); knownBox.addToRow(point1Inline); knownBox.addEndRow(point1Xline);
            knownBox.addToRow(point2X); knownBox.addToRow(point2Y); knownBox.addToRow(point2Inline); knownBox.addEndRow(point2Xline);
            knownBox.addToRow(point3X); knownBox.addToRow(point3Y); knownBox.addToRow(point3Inline); knownBox.addEndRow(point3Xline);
            add(knownBox);

            graphicBox.setBackground(Color.white);
            graphicBox.setFont(new java.awt.Font("Dialog", 1, 11));
            graphicBox.setPreferredSize(new Dimension(500, 250));
            add(graphicBox);

            buttonBox.addToRow(acceptBtn);
            buttonBox.addEndRow(cancelBtn);
            add(buttonBox);
        }

        public void accept()
        {
            if(!vols[0].checkKnownPointsOK())
            {
                new StsMessage(null, StsMessage.INFO,
                    "All fields must be specified. Fill in missing fields or cancel.");
                return;
            }
            // Fill all vols with new mapping
            if(vols.length > 1)
            {
                for(int i=1; i< vols.length; i++)
                {
                    vols[i].setPoint1Inline(vols[0].getPoint1Inline());
                    vols[i].setPoint1X(vols[0].getPoint1X());
                    vols[i].setPoint1Xline(vols[0].getPoint1Xline());
                    vols[i].setPoint1Y(vols[0].getPoint1Y());
                    vols[i].setPoint2Inline(vols[0].getPoint2Inline());
                    vols[i].setPoint2X(vols[0].getPoint2X());
                    vols[i].setPoint2Xline(vols[0].getPoint2Xline());
                    vols[i].setPoint2Y(vols[0].getPoint2Y());
                    vols[i].setPoint3Inline(vols[0].getPoint3Inline());
                    vols[i].setPoint3X(vols[0].getPoint3X());
                    vols[i].setPoint3Xline(vols[0].getPoint3Xline());
                    vols[i].setPoint3Y(vols[0].getPoint3Y());
                    vols[i].checkKnownPointsOK();
                }
            }
            parent.setVisible(false);
        }

        public void cancel()
        {
            parent.setVisible(false);
        }

        public void stateChanged(ChangeEvent e)
        {
            return;
        }
}
