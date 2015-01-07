package com.Sts.Actions.Wizards.Seismic;

import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.text.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */
public class StsViewHeadersPanel extends StsGroupBox
{
	private static final int LINE_LENGTH = 80;

	private Frame frame;
	private StsButtonFieldBean viewBinHeaderButton = new StsButtonFieldBean();
	private StsButtonFieldBean viewTxtHeaderButton = new StsButtonFieldBean();

	private StsSeismicBoundingBox[] volumes;
	private int selectedVolumeIndex;
	private DecimalFormat fmt = new DecimalFormat("#####");

	public StsViewHeadersPanel(Frame frame)
	{
		super("File Headers");
		this.frame = frame;
		constructBeans();
		jbInit();
	}

	public void constructBeans()
	{
		viewBinHeaderButton.initialize("View Binary Headers", "View the binary headers.", this, "viewBinaryHeaders");
		viewTxtHeaderButton.initialize("View Text Headers", "View the text headers.", this, "viewTextHeaders");
	}

	private void jbInit()
	{
		gbc.fill = gbc.HORIZONTAL;
		gbc.gridwidth = 3;
		addToRow(viewBinHeaderButton);
		addToRow(viewTxtHeaderButton);
	}

	public void initialize()
	{
	}

	public void setVolumes(StsSeismicBoundingBox[] volumes, int selectedVolumeIndex)
	{
		this.volumes = volumes;
		if (volumes == null || selectedVolumeIndex < 0 || selectedVolumeIndex >= volumes.length)
		{
			this.selectedVolumeIndex = -1;
		}
		else
		{
			this.selectedVolumeIndex = selectedVolumeIndex;
		}
		viewTxtHeaderButton.setEditable(! (volumes == null || volumes.length == 0));
		viewTxtHeaderButton.setEditable();
		viewBinHeaderButton.setEditable(! (volumes == null || volumes.length == 0));
		viewBinHeaderButton.setEditable();
	}

	public void viewTextHeaders()
	{
		StsTabbedTextAreaDialog dialog = new StsTabbedTextAreaDialog(frame, "SegY Text Headers", false);

		if (volumes == null)
		{
			new StsMessage(frame, StsMessage.WARNING, "No headers to display");
			return;
		}

		for (int i = 0; i < volumes.length; i++)
		{
			StsTextAreaScrollPane textArea = dialog.addTab(volumes[i].getName());
            String encoder = volumes[i].getSegyFormat().getTextHeaderFormat();
		    String header = volumes[i].segyData.readTextHdr(encoder);

			for (int j = 0; j <= header.length() - LINE_LENGTH; j += LINE_LENGTH)
			{
                String line = header.substring(j, j + LINE_LENGTH-1) + StsParameters.EOL;
				textArea.append(line);
			}
		}

		dialog.setVisible(true);
		dialog.setTabActive(selectedVolumeIndex);
	}

	public void viewBinaryHeaders()
	{
		StsTabbedTextAreaDialog dialog = new StsTabbedTextAreaDialog(frame, "SegY Binary Headers", false);
		if (volumes == null)
		{
			new StsMessage(frame, StsMessage.WARNING, "No headers to display");
			return;
		}

		for (int i = 0; i < volumes.length; i++)
		{
			StsTextAreaScrollPane textArea = dialog.addTab(volumes[i].getName());
            StsSEGYFormatRec[] allBinaryRecs = ((StsSeismicBoundingBox)volumes[i]).getSegyFormat().getAllBinaryRecords();

			for (int n = 0; n < allBinaryRecs.length; n++)
			{
				StsSEGYFormatRec rec = allBinaryRecs[n];
                textArea.append(buildLine(volumes[i].getBinaryHeaderValue(rec), rec.getDescription()));
			}
		}

		dialog.setVisible(true);
		dialog.setTabActive(selectedVolumeIndex);
	}

	private String buildLine(double num, String desc)
	{
		StringBuffer lbuf = new StringBuffer(
			"                                                                               ");

		String numText = fmt.format(num);
		if(numText.length() > 15) numText = numText.substring(0, 14);
		lbuf.insert(15 - numText.length(), numText);
		lbuf.insert(17, "-");
		lbuf.insert(19, desc);

		return lbuf.toString();
	}
}
