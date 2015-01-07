package com.Sts.Actions.Wizards.Seismic;

import com.Sts.Actions.Wizards.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.text.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.Lasseter
 * @version 1.1
 */

//TODO Not currently being used, though referenced in several places.  Should be a passive info display panel at bottom
// of wizard, but most of this info is currently being isVisible, so not sure if this is useful.
public class StsSegyInformationPanel extends StsJPanel implements ListSelectionListener
{
   private StsSEGYFormat segyFormat;
   private boolean overrideHeader = false;

   private StsGroupBox filesGroupBox = new StsGroupBox("SEGY Files to be processed");
   private StsGroupBox selectedFileBox = new StsGroupBox("Common File Properties");

   StsComboBoxFieldBean timeDepthComboBean = new StsComboBoxFieldBean();
   StsComboBoxFieldBean typeComboBean = new StsComboBoxFieldBean();
   boolean hasType = false;
   boolean tableOnly = false;
   String[] types = null;
   int currentType = 0;
   StsFloatFieldBean diskRequiredBean = new StsFloatFieldBean();
   StsIntFieldBean binaryHeaderSizeBean = new StsIntFieldBean();
   StsIntFieldBean textHeaderSizeBean = new StsIntFieldBean();
   StsButtonFieldBean viewBinHeaderButton = new StsButtonFieldBean();
   StsButtonFieldBean viewTxtHeaderButton = new StsButtonFieldBean();
   StsComboBoxFieldBean textFmtComboBean = new StsComboBoxFieldBean();
   StsFloatFieldBean startZBean = new StsFloatFieldBean();
   StsBooleanFieldBean overrideHdrBean = new StsBooleanFieldBean();
   StsIntFieldBean overrideNSamplesBean = new StsIntFieldBean();
   StsFloatFieldBean overrideSampleSpacingBean = new StsFloatFieldBean();

   private StsWizard wizard = null;

   String[] columnNames = {"Name", "Type", "Start Inline", "End Inline", "Start XLine", "End XLine",
                          "Data Min", "Data Max", "Start X", "End X",  "Start Y", "End Y"};

   private DecimalFormat numberFormat = new DecimalFormat("#,###,###,##0.0");
   private DecimalFormat fmt = new DecimalFormat("#####");
   private DecimalFormat labelFormat = new DecimalFormat("###0.0#");
   //private DecimalFormat labelFormat = new DecimalFormat("####.0#");

   StsTablePanel segyTable = new StsTablePanel(false);
   JScrollPane jScrollPane1 = new JScrollPane();

//   int selectedIndex = 0;
   byte segyFormatType = StsSEGYFormat.POSTSTACK;
   boolean editable = true;
   private StsSeismicBoundingBox[] volumes;

   private boolean giveTableVerticalHeight = true;

   JTextField messageTxt = new JTextField();

   public StsSegyInformationPanel(StsWizard wizard)
   {
      this(wizard, false, true, null);
   }

   public StsSegyInformationPanel(StsWizard wizard, boolean val)
   {
      this(wizard, val, true, null);
   }

   public StsSegyInformationPanel(StsWizard wizard, boolean val, String[] types)
   {
      this(wizard, val, true, types);
   }

   public StsSegyInformationPanel(StsWizard wizard, boolean val, String[] types, String[] cNames)
   {
      this(wizard, val, true, types);
      columnNames = cNames;
   }

   public StsSegyInformationPanel(StsWizard wizard, boolean isPrestack, boolean editable, String[] types)
   {
      if (isPrestack)
         segyFormatType = StsSEGYFormat.PRESTACK_RAW;
      if(types != null)
      {
          hasType = true;
          this.types = types;
      }
      this.editable = editable;
      this.wizard = wizard;
      constructBeans();
      jbInit();
   }

   public void setTableOnly(boolean value)
   {
       tableOnly = value;
       if(tableOnly)
           remove(selectedFileBox);
       else
           add(selectedFileBox);
       revalidate();
       repaint();
   }

   public void constructBeans()
   {
      timeDepthComboBean.initialize(wizard, "zDomainString", "Domain:", StsParameters.TD_STRINGS);
      if(hasType) typeComboBean.initialize(wizard, "objectType", "Type:", types);
      diskRequiredBean.initialize(this, "diskRequired", false, "Disk Required (MB):");
      binaryHeaderSizeBean.initialize(wizard, "binaryHeaderSize", editable, "Binary Header Size:");
      textHeaderSizeBean.initialize(wizard, "textHeaderSize", editable, "Text Header Size:");
      viewBinHeaderButton.initialize("View", "View the binary header.", this, "viewBinaryHeader");
      viewTxtHeaderButton.initialize("View", "View the text header.", this, "viewTextHeader");
      String[] textHeaderStrings = StsSEGYFormat.textHeaderFormatStrings;
      textFmtComboBean.initialize(wizard, "textHeaderFormatString", "Text Header Format:", textHeaderStrings);
      startZBean.initialize(wizard, "startZ", editable, "Start Z (+ above SL):");
      overrideHdrBean.initialize(this, "overrideHeader", "Override Header");
      overrideNSamplesBean.initialize(wizard, "overrideNSamples", false, "Samples / Trace:");
	  overrideSampleSpacingBean.initialize(wizard, "overrideSampleSpacing", false, "Sample Spacing:");
   }

   public void initialize(StsSEGYFormat segyFormat, StsWizard wizard)
   {
      this.segyFormat = segyFormat;
   }

   public void selectAllItems()
   {
       int[] rows = new int[segyTable.table.getRowCount()];
       for(int i=0; i< segyTable.table.getRowCount(); i++)
           rows[i] = i;
       segyTable.setSelectedIndices(rows);
   }

   public void clearSelectedItems()
   {
       segyTable.setSelectedIndices(null);
   }

   private void jbInit()
   {
      jScrollPane1.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
      jScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
      jScrollPane1.getViewport().add(segyTable, null);
      segyTable.addListSelectionListener(this);
      jScrollPane1.getViewport().setPreferredSize(new Dimension(560, 0));
      filesGroupBox.gbc.fill = GridBagConstraints.BOTH;
      filesGroupBox.add(jScrollPane1);

      selectedFileBox.gbc.fill = gbc.HORIZONTAL;

      StsJPanel topPanel = StsJPanel.addInsets();
      topPanel.addToRow(timeDepthComboBean);
      topPanel.addToRow(diskRequiredBean);
      if (hasType)
          topPanel.addEndRow(typeComboBean);
      selectedFileBox.add(topPanel);

      StsJPanel headerPanel = StsJPanel.addInsets();
      headerPanel.addToRow(binaryHeaderSizeBean);
      headerPanel.addToRow(viewBinHeaderButton);
      headerPanel.addEndRow(startZBean);

      headerPanel.addToRow(textHeaderSizeBean);
      headerPanel.addToRow(viewTxtHeaderButton);
      headerPanel.addEndRow(textFmtComboBean);
      selectedFileBox.add(headerPanel);

      StsJPanel overridePanel = StsJPanel.addInsets();
      overridePanel.addToRow(overrideHdrBean);
      overridePanel.addToRow(overrideNSamplesBean);
      overridePanel.addEndRow(overrideSampleSpacingBean);

      selectedFileBox.add(overridePanel);
   }

   /** As we want to control the layout of these two prefix boxes, call this method from the containing parent. */
   public void addToPanel(StsJPanel panel, boolean giveTableVerticalHeight, boolean tableOnly)
   {
       this.tableOnly = tableOnly;
       addToPanel(panel, giveTableVerticalHeight);
   }

   public void addToPanel(StsJPanel panel, boolean giveTableVerticalHeight)
   {
       panel.remove(filesGroupBox);
       panel.remove(selectedFileBox);

       panel.gbc.fill = GridBagConstraints.BOTH;
       panel.gbc.anchor = gbc.NORTH;
       if(!tableOnly)
           panel.add(selectedFileBox, 1, 0.0);
       panel.add(filesGroupBox, 1, 1.0);

       this.giveTableVerticalHeight = giveTableVerticalHeight;
   }

   public void addToPanel(StsJPanel panel)
   {
	   addToPanel(panel, true);
   }

   public boolean getOverrideHeader()
   {
      return overrideHeader;
   }

   public void viewTextHeader()
   {
      StsFileViewDialog dialog = new StsFileViewDialog(wizard.frame, "SegY Text Header View", false);
      dialog.setVisible(true);

      if (getSelectedVolume() == null)
      {
         new StsMessage(wizard.frame, StsMessage.WARNING, "Must Select PostStack3d from Table");
         return;
      }
      dialog.setViewTitle("Text Header - " + getSelectedVolume().getName());
      String encoder = segyFormat.getTextHeaderFormat();
      String header = volumes[0].segyData.readTextHdr(encoder);

      for (int i = 0; i <= header.length() - 80; i += 80)
      {
         String line = header.substring(i, i + 80);
         dialog.appendLine(line);
      }
   }

   public void viewBinaryHeader()
   {
      if (getSelectedVolume() == null)
      {
         new StsMessage(wizard.frame, StsMessage.WARNING, "Must Select PostStack3d from Table");
         return;
      }

      StsFileViewDialog dialog = new StsFileViewDialog(wizard.frame, "SegY Binary Header View", false);
      dialog.setVisible(true);
      dialog.setViewTitle("Binary Header - " + getSelectedVolume().getName());

      StsSeismicBoundingBox selectedVolume = getSelectedVolume();

      StsSEGYFormatRec[] allBinaryRecs = null;
      allBinaryRecs = segyFormat.getAllBinaryRecords();

      for (int n = 0; n < allBinaryRecs.length; n++)
      {
         StsSEGYFormatRec rec = allBinaryRecs[n];
         dialog.appendLine(buildLine(selectedVolume.getBinaryHeaderValue(rec), rec.getDescription()));
      }
   }

   public void setSegyFormat(StsSEGYFormat format)
   {
       segyFormat = format;
   }

   public void valueChanged(ListSelectionEvent e)
   {
       if(!e.getValueIsAdjusting())
       {
           wizard.itemSelected();
           setSelectedValues();
       }
   }

   public void setSelectedValues()
   {
       binaryHeaderSizeBean.getValueFromPanelObject();
       textHeaderSizeBean.getValueFromPanelObject();
       startZBean.getValueFromPanelObject();
       overrideNSamplesBean.getValueFromPanelObject();
       overrideSampleSpacingBean.getValueFromPanelObject();
       timeDepthComboBean.getValueFromPanelObject();

       if(hasType)
          typeComboBean.getValueFromPanelObject();

       textFmtComboBean.getValueFromPanelObject();
       overrideHdrBean.getValueFromPanelObject();
   }

   public String buildLine(double num, String desc)
   {
      StringBuffer lbuf = new StringBuffer(
         "                                                                               ");

      String numText = fmt.format(num);
      lbuf.insert(10 - numText.length(), numText);
      lbuf.insert(12, "-");
      lbuf.insert(14, desc);

      return lbuf.toString();
   }

   public void setTableValues(StsSeismicBoundingBox[] volumes_)
   {
      final StsSeismicBoundingBox[] volumes = volumes_;
      Runnable runnable = new Runnable()
      {
        public void run()
        {
            runSetTableValues(volumes);
        }
      };
      StsToolkit.runWaitOnEventThread(runnable);
   }

   private void runSetTableValues(StsSeismicBoundingBox[] volumes)
   {
      // Need to allocate and fill based on user specified attributes - SAJ
      this.volumes = volumes;
      int[] selectedRows = segyTable.getSelectedIndices();
      segyTable.removeListSelectionListener(this);
	  segyTable.setTableModel(new InformationTableModel(volumes, columnNames));
	  int height = 0;
	  if (volumes != null)
	  {
		  height = volumes.length * 25;
          if(height > 100)
              height = 100;
	  }
      segyTable.setPreferredSize(new Dimension(11 * 75, height));
      segyTable.setSelectedIndices(selectedRows);
      segyTable.addListSelectionListener(this);
	  if(giveTableVerticalHeight && volumes != null)
	  {
		  jScrollPane1.getViewport().setPreferredSize(new Dimension(560, height));
	  }
	  else
	  {
		  jScrollPane1.getViewport().setPreferredSize(new Dimension(560, 0));
	  }
   }

   public int getTypeFromString(String value)
   {
       for(int i=0; i<types.length; i++)
           if(value.equals(types[i]))
               return i;
       return 0;
   }

   public String getStringFromType(int type)
   {
       if(hasType)
           return types[type];
       else
           return "";
   }

   public float getDiskRequired()
   {
      return StsSeismicBoundingBox.calcDiskRequired(volumes);
   }

   public void setOverrideNSamples(int value)
   {
	   // this will be called in the wizard
      //segyFormat.setOverrideNSamples(value);
      overrideNSamplesBean.setValue(value);
   }

   public void setOverrideSampleSpacing(float value)
   {
      overrideSampleSpacingBean.setValue(value);
   }

   public void setEditable(boolean value)
   {
      timeDepthComboBean.setEditable(value);
   }

   public void setSelectedIndex(int idx)
   {
      segyTable.setSelectedIndices(new int[] {idx});
      System.out.println("setSelectedIndex:" + idx);

      return;
   }

   public void setSelectedIndices(int[] indices)
   {
       segyTable.setSelectedIndices(null);
       if(indices.length == 0)
           return;
       segyTable.setSelectedIndices(indices);
       System.out.println("setSelectedIndices:" + indices[0] + ":" + indices[indices.length-1]);

       return;
   }

   /*
       public void setSelectedVolume(StsRotatedGridBoundingBox vol)
       {
           for(int i=0; i<segyTable.getNumberOfRows(); i++)
           {
               if(vol.getName().equals(volumes[i].getName()))
                  segyTable.setSelectedIndices(new int[] {i});
           }
           setValues(vol);
           return;
       }
    */
   public StsSeismicBoundingBox getSelectedVolume()
   {
      if (volumes == null) return null;
      int selectedIndex = 0;
      int[] selectedIndices = segyTable.getSelectedIndices();
      if (selectedIndices != null && selectedIndices.length > 0)
         selectedIndex = selectedIndices[0];
      if (volumes.length < selectedIndex + 1)
         return null;
      return volumes[selectedIndex];
   }

   public StsSeismicBoundingBox[] getSelectedVolumes()
   {
      if (volumes == null) return null;
      int[] selectedIndices = segyTable.getSelectedIndices();
      if(selectedIndices.length == 0)
          return null;

      StsSeismicBoundingBox[] selectedVolumes = new StsSeismicBoundingBox[selectedIndices.length];
      for(int i=0; i<selectedVolumes.length; i++)
      {
          selectedVolumes[i] = volumes[selectedIndices[i]];
      }
      return selectedVolumes;
   }

   public void setMessage(String msg)
   {
      messageTxt.setText(msg);
      repaint();
   }

   static public void main(String[] args)
   {
      com.Sts.MVC.StsModel model = new com.Sts.MVC.StsModel();
      com.Sts.MVC.StsProject project = new com.Sts.MVC.StsProject();
      model.setProject(project);
      StsSegyInformationPanel panel = new StsSegyInformationPanel(null, true, true, null);
      StsSEGYFormat segyFormat = StsSEGYFormat.constructor(model, StsSEGYFormat.PRESTACK_RAW);
      panel.initialize(segyFormat, null);
      StsJPanel backPanel = StsJPanel.addInsets();
      panel.addToPanel(backPanel, true);
      StsToolkit.createDialog(backPanel, true);
   }

   private class InformationTableModel extends AbstractTableModel
   {
	   private StsSeismicBoundingBox[] volumes = null;
	   String[] columnNames = null;

       public InformationTableModel(StsSeismicBoundingBox[] volumes, String[] cNames)
       {
           super();
           this.volumes = volumes;
           this.columnNames = cNames;
       }

	   public int getColumnCount()
	   {
		   if (volumes == null)
		   {
			   return 0;
		   }
		   return columnNames.length;
	   }

	   public String getColumnName(int columnIndex)
	   {
		   return columnNames[columnIndex];
	   }

	   public boolean isCellEditable(int rowIndex, int columnIndex)
	   {
		   return false;
	   }

	   public Class getColumnClass(int columnIndex)
	   {
		   return String.class;
	   }

	   public int getRowCount()
	   {
		   if (volumes == null)
		   {
			   return 0;
		   }
		   return volumes.length;
	   }

	   public Object getValueAt(int rowIndex, int columnIndex)
	   {
		   switch (columnIndex)
		   {
			   case 0:
				   return volumes[rowIndex].getName();
               case 1:
				   return volumes[rowIndex].getTypeAsString();
			   case 2:
				   return labelFormat.format(new Float(volumes[rowIndex].getRowNumMin()));
			   case 3:
				   return labelFormat.format(new Float(volumes[rowIndex].getRowNumMax()));
			   case 4:
				   return labelFormat.format(new Float(volumes[rowIndex].getColNumMin()));
			   case 5:
				   return labelFormat.format(new Float(volumes[rowIndex].getColNumMax()));
			   case 6:
				   return labelFormat.format(new Float(volumes[rowIndex].getDataMin()));
			   case 7:
				   return labelFormat.format(new Float(volumes[rowIndex].getDataMax()));
			   case 8:
				   return labelFormat.format(volumes[rowIndex].getXOrigin());
			   case 9:
				   return labelFormat.format(volumes[rowIndex].getXOrigin() + (volumes[rowIndex].getXInc() * volumes[rowIndex].getXSize()));
			   case 10:
				   return labelFormat.format(volumes[rowIndex].getYOrigin());
			   case 11:
				   return labelFormat.format(volumes[rowIndex].getYOrigin() + (volumes[rowIndex].getYInc() * volumes[rowIndex].getYSize()));
		   }
		   return "----";
	   }
   }
}
