package com.Sts.Actions.Wizards.Seismic;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.Lasseter
 * @version 1.1
 */

public class StsSeismic3dRangeEditPanel extends StsFieldBeanPanel
{
    public CroppedBoundingBox croppedBoundingBox;
    /** the volumes selected for this pass thru the wizard */
    StsSeismicBoundingBox[] selectedVolumes;

    private boolean rowNumInverted = false;

    private StsGroupBox rangeBox = new StsGroupBox("Editable Post-Stack Information");
    private StsGroupBox originAngleBox = new StsGroupBox("Other Post-Stack Information");
    private StsFloatFieldBean rowNumMinBean;
    private StsFloatFieldBean rowNumMaxBean;
    private StsFloatFieldBean rowNumIncBean;
    private StsFloatFieldBean colNumMinBean;
    private StsFloatFieldBean colNumMaxBean;
    private StsFloatFieldBean colNumIncBean;
    private StsFloatFieldBean zMinBean;
    private StsFloatFieldBean zMaxBean;
    private StsFloatFieldBean zIncBean;
	private JLabel minLabel = new JLabel("Minimum");
	private JLabel maxLabel = new JLabel("Maximum");
	private JLabel intervalLabel = new JLabel("Increment");
    private StsDoubleFieldBean xOriginBean;
    private StsDoubleFieldBean yOriginBean;
    private StsFloatFieldBean angleBean;


    public StsSeismic3dRangeEditPanel()
	{
        try
		{
            croppedBoundingBox = new CroppedBoundingBox();
            constructBeans();
            buildPanel();
        }
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

	public void initialize(StsSeismicBoundingBox[] volumes)
	{
        croppedBoundingBox.initialize(volumes);
        croppedBoundingBox.initializeBeanValuesAndRanges();
        setPanelObject(croppedBoundingBox);
    }

    private void constructBeans()
    {
        rowNumMinBean = new StsFloatFieldBean(CroppedBoundingBox.class, "rowNumMin", "Inline:");
        rowNumMaxBean = new StsFloatFieldBean(CroppedBoundingBox.class, "rowNumMax");
        rowNumIncBean = new StsFloatFieldBean(CroppedBoundingBox.class, "rowNumInc");
        colNumMinBean = new StsFloatFieldBean(CroppedBoundingBox.class, "colNumMin", "Crossline:");
        colNumMaxBean = new StsFloatFieldBean(CroppedBoundingBox.class, "colNumMax");
        colNumIncBean = new StsFloatFieldBean(CroppedBoundingBox.class, "colNumInc");
        zMinBean      = new StsFloatFieldBean(CroppedBoundingBox.class, "zMin", "Time/Depth:");
        zMaxBean      = new StsFloatFieldBean(CroppedBoundingBox.class, "zMax");
        zIncBean      = new StsFloatFieldBean(CroppedBoundingBox.class, "zInc");
        xOriginBean   = new StsDoubleFieldBean(CroppedBoundingBox.class, "xOrigin", false, "X Origin:");
        yOriginBean   = new StsDoubleFieldBean(CroppedBoundingBox.class, "yOrigin", false, "Y Origin:");
        angleBean     = new StsFloatFieldBean(CroppedBoundingBox.class, "angle", false, "Angle:");
    }

    private void buildPanel()
	{
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.NORTH;

        rangeBox.gbc.gridx = 1; // leave empty cell
        rangeBox.addToRow(minLabel);
        rangeBox.addToRow(maxLabel);
        rangeBox.addEndRow(intervalLabel);

        rangeBox.addToRow(rowNumMinBean);
        rangeBox.addToRow(rowNumMaxBean);
        rangeBox.addEndRow(rowNumIncBean);

        rangeBox.addToRow(colNumMinBean);
        rangeBox.addToRow(colNumMaxBean);
        rangeBox.addEndRow(colNumIncBean);

        rangeBox.addToRow(zMinBean);
        rangeBox.addToRow(zMaxBean);
        rangeBox.addEndRow(zIncBean);
        
        add(rangeBox);

        originAngleBox.addToRow(xOriginBean);
        originAngleBox.addToRow(yOriginBean);
        originAngleBox.addEndRow(angleBean);
        add(originAngleBox);
    }

    public void disableCropOptions()
    {
        rowNumMinBean.setEditable( false);
        rowNumMaxBean.setEditable( false);
        rowNumIncBean.setEditable( false);
        colNumMinBean.setEditable( false);
        colNumMaxBean.setEditable( false);
        colNumIncBean.setEditable( false);
        zMinBean.setEditable( false);
        zMaxBean.setEditable( false);
        zIncBean.setEditable( false);
    }

    public void initializeBeanValuesAndRanges()
    {
        croppedBoundingBox.initializeBeanValuesAndRanges();
    }

    class CroppedBoundingBox extends StsCroppedBoundingBox
    {
        CroppedBoundingBox()
        {
            super(false);
        }

        private void initialize(StsSeismicBoundingBox[] volumes)
        {
            StsSeismicBoundingBox firstVolume = volumes[0];
            initializeOriginAndAngle(firstVolume);

            for(int n = 0; n < volumes.length; n++)
                addBoundingBox(volumes[n]);

            nSlices = volumes[0].nSlices;
            checkCurrentRangeValues();
            originalBoundingBox = new StsRotatedGridBoundingBox(this);
        }

        /** called when rangeEditPanel is initialized, or one of the range values is changed
         *  in order to check that ranges are ok; if so enable nextButton; otherwise disable it.
         */
        private boolean checkCurrentRangeValues()
        {
            //todo should do a more thorough range check
            return nSlices > 0 && zInc > 0.0f;
        }

        private void initializeBeanValuesAndRanges()
        {
            rowNumInverted = rowNumMin > rowNumMax;
            if(!rowNumInverted)
            {
                rowNumMinBean.setValueAndRangeFixStep(rowNumMin, rowNumMin, rowNumMax, rowNumInc);
                rowNumMaxBean.setValueAndRangeFixStep(rowNumMax, rowNumMin, rowNumMax, rowNumInc);
                rowNumIncBean.setValueAndRangeFixStep(rowNumInc, rowNumInc, rowNumMax - rowNumMin, rowNumInc);

            }
            else
            {
                rowNumMinBean.setValueAndRangeFixStep(rowNumMax, rowNumMax, rowNumMin, -rowNumInc);
                rowNumMaxBean.setValueAndRangeFixStep(rowNumMin, rowNumMax, rowNumMin, -rowNumInc);
                rowNumIncBean.setValueAndRangeFixStep(-rowNumInc, -rowNumInc, rowNumMin - rowNumMax, -rowNumInc);
            }

            colNumMinBean.setValueAndRangeFixStep(colNumMin, colNumMin, colNumMax, colNumInc);
            colNumMaxBean.setValueAndRangeFixStep(colNumMax, colNumMin, colNumMax, colNumInc);
            colNumIncBean.setValueAndRangeFixStep(colNumInc, colNumInc, colNumMax - colNumMin, colNumInc);
        }

        public void setRowNumMin(float rowNumMin)
        {
            if(!rowNumInverted)
                this.rowNumMin = rowNumMin;
            else
                this.rowNumMax = rowNumMin;
        }
        public float getRowNumMin()
        {
            if(!rowNumInverted)
                return rowNumMin;
            else
                return rowNumMax;
        }

        public void setRowNumMax(float rowNumMax)
        {
            if(!rowNumInverted)
                this.rowNumMax = rowNumMax;
            else
                this.rowNumMin = rowNumMax;
        }

        public float getRowNumMax()
        {
            if(!rowNumInverted)
                return rowNumMax;
            else
                return rowNumMin;
        }

        public void setRowNumInc(float value)
        {
//            rowNumMinBean.fixStep(value);
//            rowNumMaxBean.fixStep(value);
        	if(rowNumInc == value)	return;
            rowNumInc = value;
            if(!rowNumInverted)
            {
                rowNumMax = StsMath.floor((originalBoundingBox.rowNumMax - rowNumMin) / rowNumInc) * rowNumInc + rowNumMin;
                rowNumMaxBean.setValue(rowNumMax);
            }
            else
            {
                rowNumMin = StsMath.floor((originalBoundingBox.rowNumMin - rowNumMax) / rowNumInc) * rowNumInc + rowNumMax;
                rowNumMaxBean.setValue(rowNumMin);
                rowNumInc = -value;
            }
        }

        public float getRowNumInc()
        {
            if(!rowNumInverted)
                return rowNumInc;
            else
                return -rowNumInc;
        }

        public void setColNumInc(float value)
        {
        	if(colNumInc == value) return;
            colNumInc = value;
//            colNumMinBean.fixStep(colNumInc);
//            colNumMaxBean.fixStep(colNumInc);
            colNumMax = StsMath.floor((colNumMax - colNumMin) / colNumInc) * colNumInc + colNumMin;
            colNumMaxBean.setValue(colNumMax);
        }

        public void setZMin(float value)
        {
            zMin = value;
            zMax = zMin + (nSlices-1)*zInc;
            zMax = Math.min(zMax, originalBoundingBox.zMax);
            zMaxBean.setValue(zMax);
            checkCurrentRangeValues();
        }

        public void setZMax(float value)
        {
            zMax = value;
//            zMaxBean.setValue(zMax);
//            zMin = zMax - (nCroppedSlices-1)*zInc;
//            zMinBean.setValue(zMin);
        }

        public void setZInc(float value)
        {
            zInc = value;
            zMinBean.fixStep(zInc);
            zMaxBean.fixStep(zInc);
            if(zMax > zMin)
                nSlices = (int)((zMax - zMin)/zInc) + 1;
            zMax = zMin + (nSlices-1)*zInc;
            zMaxBean.setValue(zMax);
            checkCurrentRangeValues();
        }
    }

    static public void main(String[] args)
	{
		StsSeismic3dRangeEditPanel panel = new StsSeismic3dRangeEditPanel();
        StsToolkit.createDialog(panel);
    }
}