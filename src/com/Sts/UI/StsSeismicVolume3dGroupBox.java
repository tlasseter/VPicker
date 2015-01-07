package com.Sts.UI;

import com.Sts.DBTypes.StsSeismicVelocityModel;
import com.Sts.DBTypes.StsSeismicVolume;
import com.Sts.MVC.StsModel;
import com.Sts.Types.StsSeismicBoundingBox;
import com.Sts.UI.Beans.StsFloatFieldBean;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Beans.StsStringFieldBean;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsParameters;
import com.Sts.Utilities.StsToolkit;

import javax.swing.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: 3/2/11
 * Time: 9:55 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsSeismicVolume3dGroupBox extends StsSeismicVolume2dGroupBox
{
    StsFloatFieldBean rowNumMinBean;
    StsFloatFieldBean rowNumMaxBean;
    StsFloatFieldBean rowNumIncBean;
    StsFloatFieldBean colNumMinBean;
    StsFloatFieldBean colNumMaxBean;
    StsFloatFieldBean colNumIncBean;

    public StsSeismicVolume3dGroupBox(StsSeismicBoundingBox volume, String title)
    {
  		super(volume, title);
    }

    public StsSeismicVolume3dGroupBox(StsSeismicBoundingBox volume, String title, StsSeismicBoundingBox originalBoundingBox)
    {
  		super(volume, title, originalBoundingBox);
    }

    protected void buildBeans(boolean editable)
    {
        super.buildBeans(editable);
		rowNumMinBean = new StsFloatFieldBean(volume, "rowNumMin", editable);
		rowNumMaxBean = new StsFloatFieldBean(volume, "rowNumMax", editable);

		rowNumIncBean = new StsFloatFieldBean(volume, "rowNumInc",editable, null, true);
        float rowNumInc = volume.rowNumInc;
        rowNumIncBean.setRangeFixStep(rowNumInc, 100*rowNumInc, rowNumInc);

		colNumMinBean = new StsFloatFieldBean(volume, "colNumMin", editable);
		colNumMaxBean = new StsFloatFieldBean(volume, "colNumMax", editable);
		colNumIncBean = new StsFloatFieldBean(volume, "colNumInc", editable, null, true);
		float colNumInc = volume.colNumInc;
        colNumIncBean.fixStep(volume.colNumInc);
        if(colNumInc > 0.0f)
            colNumIncBean.setRange(volume.colNumInc, 100*volume.colNumInc);
        else
            colNumIncBean.setRange(100*volume.colNumInc, volume.colNumInc);
    }

    protected void addCoordinatesPanel()
    {
        StsJPanel coordinatesPanel = new StsJPanel();

        coordinatesPanel.gbc.anchor = WEST;

        coordinatesPanel.addToRow(new JLabel("Inline:"));
		coordinatesPanel.addToRow(rowNumMinBean);
		coordinatesPanel.addToRow(rowNumMaxBean);
		coordinatesPanel.addEndRow(rowNumIncBean);

        coordinatesPanel.addToRow(new JLabel("Crossline:"));
		coordinatesPanel.addToRow(colNumMinBean);
		coordinatesPanel.addToRow(colNumMaxBean);
		coordinatesPanel.addEndRow(colNumIncBean);

        coordinatesPanel.addToRow(zOrTLabel);
		coordinatesPanel.addToRow(zMinBean);
		coordinatesPanel.addToRow(zMaxBean);
		coordinatesPanel.addEndRow(zIncBean);

        add(coordinatesPanel);
    }

    public void adjustZRange(String zDomainString)
    {
        if(zDomainString  == StsParameters.TD_DEPTH_STRING) // originalVolume.zDomain == TD_TIME
		{
			float zMin = originalVolume.getMinDepthAtTime(originalVolume.zMin);
			float zMax = originalVolume.getMaxDepthAtTime(originalVolume.zMax);
			int nVolumeSlices = originalVolume.getNSlices();
			double[] niceScale = StsMath.niceScale(zMin, zMax, nVolumeSlices, true);
			zMin = (float)niceScale[0];
			zMax = (float)niceScale[1];
			float zInc = (float)niceScale[2];
			nVolumeSlices = Math.round((zMax - zMin)/zInc);

            volume.zMin = zMin;
            volume.zMax = zMax;
            volume.zInc = zInc;
            volume.nSlices = nVolumeSlices;

            zMinBean.setValue(zMin);
		    zMaxBean.setValue(zMax);
		    zIncBean.setValue(zInc);
		}
        //TODO need to implement getMinTimeAtDepth and getMaxTimeAtDepth
        else // exportDomainString == StsParameters.TD_TIME_STRING) && volumeDomainString == TD_DEPTH_STRING
        {

        }
	}

    static public void main(String[] args)
    {
        try
        {
            StsSeismicVolume volume = StsSeismicVolume.constructTestVolume(args[0]);
            StsSeismicVolume3dGroupBox groupBox = new StsSeismicVolume3dGroupBox(volume, "Test Seismic Volume Group Box");
            StsToolkit.createDialog("Test", groupBox, true, 250, 500);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsSeismicVolume3dGroupBox.class, "main", e);
        }
    }
}
