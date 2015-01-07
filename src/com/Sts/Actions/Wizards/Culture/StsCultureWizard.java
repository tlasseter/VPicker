package com.Sts.Actions.Wizards.Culture;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import org.kabeja.dxf.*;
import org.kabeja.parser.*;

import java.awt.*;
import java.io.*;
import java.util.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class StsCultureWizard extends StsWizard
{
    double pi=Math.PI/180;
    int bx=0,by=0;
    int nvert=0,maxvert=100;
    float THX=0,THY=10;
    float vert[][]=new float[maxvert][10],light[] = {1,1,1};
    String str=null; // filename
    Image img=null;
    boolean bo=false;
    Scrollbar sb=null;
    
    StsAbstractFileSet fileSet;
    TreeMap filenameTreeMap;
    ObjectFilenameFilter filenameFilter;
    StsCultureObjectSet2D cultureSet2D = null;
    String cultureSetName = new String("Culture");

    public StsCultureSetSelect cultureSetSelect = new StsCultureSetSelect(this);
    public StsCultureSetAttributes cultureSetAttributes = new StsCultureSetAttributes(this);
    public StsCultureSetLoad cultureSetLoad = new StsCultureSetLoad(this);
    private StsWizardStep[] mySteps =
    {
        cultureSetSelect, cultureSetAttributes, cultureSetLoad
    };

    public StsCultureWizard(StsActionManager actionManager)
    {
        super(actionManager,500,500);
        addSteps(mySteps);
        filenameFilter = new ObjectFilenameFilter(new String[] {"xml","dxf"});
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Process & Load Culture");
        String dirname = new String(model.getProject().getRootDirString() + "/");
        initializeCultureFilenameSet(dirname);
        disableFinish();
        return super.start();
    }

    public boolean end()
    {
        if (success) model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == cultureSetAttributes)
            cultureSetLoad.constructPanel();
        gotoNextStep();
    }

    public void finish()
    {
        super.finish();
    }

    public boolean initializeCultureFilenameSet(String dirname)
    {
        try
        {
            fileSet = StsFileSet.constructor(dirname, filenameFilter);
            filenameTreeMap = new TreeMap();
            StsAbstractFile[] files = fileSet.getFiles();
            for (int n = 0; n < files.length; n++)
            {
                String filename = files[n].getFilename();
                String objectName = filenameFilter.getFilenameStem(filename);
                filenameTreeMap.put(objectName, files[n]);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsCultureWizard.getCultureFilenameSets() failed.",
                                         e, StsException.WARNING);
            return false;
        }
    }

    public Object[] getObjectNames()
    {
        return filenameTreeMap.keySet().toArray();
    }

    public StsAbstractFile[] getSelectedFiles()
    {
        Object[] objectNames = cultureSetSelect.getSelectedObjects();
        int nSelected = objectNames.length;
        StsAbstractFile[] selectedFiles = new StsAbstractFile[nSelected];
        for (int n = 0; n < nSelected; n++)
            selectedFiles[n] = (StsAbstractFile)filenameTreeMap.get(objectNames[n]);
        return selectedFiles;
    }

    final class ObjectFilenameFilter implements FilenameFilter
    {
        String[] filters = null;
        int length;

        public ObjectFilenameFilter(String[] filters)
        {
        	this.filters = new String[filters.length];
        	for(int i=0; i<filters.length; i++)
        		this.filters[i] = filters[i].toLowerCase();
        }

        public boolean accept(File dir, String name)
        {
        	for(int i=0; i<filters.length; i++)
        	{
        		if(name.toLowerCase().endsWith(filters[i]))
        			return true;
        	}
        	return false;
        }

        public String getFilenameStem(String filename)
        {
        	for(int i=0; i<filters.length; i++)
        	{        	
        		int filterStart = filename.toLowerCase().indexOf(filters[i]);
        		if(filterStart > 0)
        			return filename.substring(0, filterStart - 1);
        	}
            return null;
        }
    }

    public boolean createCultureSet(StsProgressPanel panel)
    {
       if( model == null ) return false;

       try
       {
           StsAbstractFile[] selectedFiles = getSelectedFiles();
           model.disableDisplay();

           disablePrevious();

           // turn on the wait cursor
           StsCursor cursor = new StsCursor(model.win3d, Cursor.WAIT_CURSOR);
           int nSelected = selectedFiles.length;
           int nLoaded = 0;

           cultureSet2D = new StsCultureObjectSet2D(false);
           cultureSet2D.setName(getCultureSetName());
           StsProject project = model.getProject();
           cultureSet2D.addCultureDisplayable(project);
           panel.initialize(nSelected);
           String fileType = filenameFilter.filters[0];
           for(int n = 0; n < nSelected; n++)
           {
               if (selectedFiles[n] == null) continue;
               
               for(int j=0; j<filenameFilter.filters.length; j++)
               {
            	   if(selectedFiles[n].getFilename().toLowerCase().endsWith(filenameFilter.filters[j]))
            		  fileType = filenameFilter.filters[j];
               }
               panel.appendLine("Loading " + fileType + " file named: " + selectedFiles[n].getFilename());
               if(!loadCultureObject(selectedFiles[n], panel, fileType))
               {
            	   panel.appendLine("Failed to load file named: " + selectedFiles[n].getFilename());
            	   panel.setLevel(StsProgressBar.WARNING);
               }
               else
            	   nLoaded++;
               panel.setValue(n+1);
               panel.setDescription("Culture set " + nLoaded + " of " + nSelected + " loaded.");
           }
           panel.setDescription("Loaded " + nLoaded + " culture files of " + nSelected + " selected.");
           cultureSet2D.setCultureDisplayable(model.getProject());
           
           // turn off the wait cursor
           if( cursor!=null ) cursor.restoreCursor();

           disableCancel();
           enableFinish();
           
           if(cultureSet2D.getNumberOfObjects() > 0)
           {
        	   cultureSet2D.addToProject();
        	   project.adjustBoundingBoxes(true, false);
        	   project.checkAddUnrotatedClass(StsCultureObjectSet2D.class);
        	   project.rangeChanged(); 
        	   cultureSet2D.addToModel();           
           
        	   model.win3d.cursor3d.initialize();
        	   model.win3d.cursor3dPanel.setSliderValues();
        	   model.getGlPanel3d().setDefaultView();          
        	   model.enableDisplay();
        	   model.win3dDisplay(); // display the surfaces
           }

           return true;
       }
       catch (Exception e)
       {
           panel.appendLine("StsCultureWizard.createCultureSets() failed." + e.getMessage());
           panel.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
           StsException.outputException("StsCultureWizard.createCultureSets() failed.", e, StsException.WARNING);
           return false;
       }
    }

    private boolean loadCultureObject(StsAbstractFile file, StsProgressPanel panel, String fileType)
    {
        if(fileType.equals("dxf"))
        	return loadDxfCulture(file, panel);
        else
        	return loadXmlCulture(file, panel);
    }
    
    private boolean loadDxfCulture(StsAbstractFile file, StsProgressPanel panel)
    {
        StsCultureLine2D line = null; 
        StsCulturePointSet2D pointSet = null;
        StsCultureText text = null;
        
        byte type = -1;
        int size = 1;
        int stroke = 0;
        int nLoaded = 0;
        
        Parser parser = ParserBuilder.createDefaultParser();
        try 
        {
            InputStream in = null;
            try
            {
                in = new FileInputStream(file.getPathname()); // Does not handle spaces???
            } 
            catch (FileNotFoundException e)
            {
                System.out.println("File not found....");
                return false;
            }

            StsColor color;
            StsSpectrum autoCadSpectrum = model.getSpectrum(StsSpectrumClass.SPECTRUM_AUTOCAD);
            
            StsColor overrideColor = cultureSetAttributes.panel.getStsColor();
            boolean override = cultureSetAttributes.panel.getOverrideColor();
            float coordinateScalar = cultureSetAttributes.panel.getScalar();

            //parse
            parser.parse(in, DXFParser.DEFAULT_ENCODING);
            //get the documnet and the layer
            DXFDocument doc = parser.getDocument();
            // iterate over layers
            Iterator iter = doc.getDXFLayerIterator();
            while(iter.hasNext())
            {
            	DXFLayer layer = (DXFLayer)iter.next();
            	color = autoCadSpectrum.getColor(layer.getColor());
            	if(override) color = overrideColor;
            	java.util.List face3dList = layer.getDXFEntities(DXFConstants.ENTITY_TYPE_3DFACE);
            	if(face3dList != null)
            	{
            		Iterator iter2 = face3dList.iterator();
            		while(iter2.hasNext())
            		{
            			DXF3DFace face3d = (DXF3DFace)iter2.next();
            			color = autoCadSpectrum.getColor(face3d.getColor());
                    	if(override) color = overrideColor;
            			line = new StsCultureLine2D();     
            			line.setType(StsCultureObject2D.LINE);
            			line.setSize(1);
            			line.setStsColor(color); // Color of Line
            			line.setStroke(StsCultureLine2D.SOLID);
                    
            			line.addPoint(face3d.getPoint1().getX() * coordinateScalar, face3d.getPoint1().getY() * coordinateScalar, face3d.getPoint1().getZ(), model.getProject().getTimeMin());
            			line.addPoint(face3d.getPoint2().getX() * coordinateScalar, face3d.getPoint2().getY() * coordinateScalar, face3d.getPoint2().getZ(), model.getProject().getTimeMin());
            			line.addPoint(face3d.getPoint3().getX() * coordinateScalar, face3d.getPoint3().getY() * coordinateScalar, face3d.getPoint3().getZ(), model.getProject().getTimeMin());
            			line.addPoint(face3d.getPoint4().getX() * coordinateScalar, face3d.getPoint4().getY() * coordinateScalar, face3d.getPoint4().getZ(), model.getProject().getTimeMin());
            			cultureSet2D.addObject(line); 
            			nLoaded++;
            		}
            	}
            	java.util.List polylineList = layer.getDXFEntities(DXFConstants.ENTITY_TYPE_POLYLINE);
            	if(polylineList != null)
            	{
            		Iterator iter2 = polylineList.iterator();
            		while(iter2.hasNext())
            		{
            			DXFPolyline polyline = (DXFPolyline)iter2.next();
            			color = autoCadSpectrum.getColor(polyline.getColor());
                    	if(override) color = overrideColor;
            			line = new StsCultureLine2D();     
            			line.setType(StsCultureObject2D.LINE);
            			line.setSize(1);
            			line.setStsColor(color); // Color of Line
            			line.setStroke(StsCultureLine2D.SOLID);
                    
            			Iterator iter3 = polyline.getVertexIterator();
            			while(iter3.hasNext())
            			{
            				DXFVertex vertex = (DXFVertex)iter3.next();
            				line.addPoint(vertex.getX() * coordinateScalar, vertex.getY() * coordinateScalar, vertex.getZ(), model.getProject().getTimeMin());
            			}
            			cultureSet2D.addObject(line);       
            			nLoaded++;
            		}
            	}
            	java.util.List pointList = layer.getDXFEntities(DXFConstants.ENTITY_TYPE_POINT);           	
            	if(pointList != null)
            	{
            		Iterator iter2 = pointList.iterator();
        			pointSet = new StsCulturePointSet2D();
        			pointSet.setType(StsCulturePointSet2D.POINT);
        			
            		while(iter2.hasNext())
            		{
            			DXFPoint point = (DXFPoint)iter2.next();
            			color = autoCadSpectrum.getColor(point.getColor());
                    	if(override) color = overrideColor;
            			pointSet.addPoint(point.getX() * coordinateScalar,point.getY() * coordinateScalar,point.getZ(),model.getProject().getTimeMin());
            			cultureSet2D.addObject(pointSet);
            			nLoaded++;
            		}
            		pointSet.setStsColor(color);
            	}            	
            }
            
            // iterate over blocks
            iter = doc.getDXFBlockIterator();
            DXFBlock block = null;
            int count = 1;
            while(iter.hasNext())
            {
            	block = (DXFBlock)iter.next();
            	Iterator iter2 = block.getDXFEntitiesIterator();
            	while(iter2.hasNext())
            	{
            		Object dxfEntity = iter2.next();
            		if(dxfEntity instanceof DXFPolyline)
            		{
            			color = autoCadSpectrum.getColor(((DXFPolyline)dxfEntity).getColor());
                    	if(override) color = overrideColor;
                        line = new StsCultureLine2D();     
                        line.setType(StsCultureObject2D.LINE);
                        line.setSize(1);
                        line.setStsColor(color); // Color of Line
                        line.setStroke(StsCultureLine2D.SOLID);
                        
            			Iterator iter3 = ((DXFPolyline)dxfEntity).getVertexIterator();
            			while(iter3.hasNext())
            			{
            				DXFVertex vertex = (DXFVertex)iter3.next();
            				line.addPoint(vertex.getX() * coordinateScalar, vertex.getY() * coordinateScalar, vertex.getZ(), model.getProject().getTimeMin());
            			}
            			cultureSet2D.addObject(line);
            			nLoaded++;
            		}           		
            		else if(dxfEntity instanceof DXFText)
            		{
            			DXFText dxfText = (DXFText)dxfEntity;
            			color = autoCadSpectrum.getColor(dxfText.getColor());
                    	if(override) color = overrideColor;
                        text = new StsCultureText();
                        type = StsCultureObject2D.TEXT;
                        text.setSize(size);
                        text.setStsColor(color); // Color of Point;
                        text.setText(dxfText.getText());
                        text.setPoint(dxfText.getAlignX() * coordinateScalar, dxfText.getAlignY() * coordinateScalar, dxfText.getAlignZ(), model.getProject().getTimeMin());
                        cultureSet2D.addObject(text);
                        nLoaded++;
            		}
            		else if(dxfEntity instanceof DXFPoint)
            		{
            			DXFPoint point = (DXFPoint)dxfEntity;
            			color = autoCadSpectrum.getColor(point.getColor());
                    	if(override) color = overrideColor;
            			pointSet = new StsCulturePointSet2D();
            			pointSet.setType(StsCulturePointSet2D.POINT);
            			
                		pointSet.addPoint(point.getX() * coordinateScalar,point.getY() * coordinateScalar,point.getZ(),model.getProject().getTimeMin());
                		cultureSet2D.addObject(pointSet); 
                		nLoaded++;
            		}
            		else
            			StsMessageFiles.infoMessage("Unsupportted DXF Entity type (" + dxfEntity.toString() + ")");
            	}
            }
        } 
        catch (Exception e) 
        {
        	StsException.outputException("Unable to load selected DXF file", e, StsException.FATAL);
        	e.printStackTrace();
        	return false;
        }
        if(nLoaded == 0) 
        	return false;
        else
        	return true;
    }

    private boolean loadXmlCulture(StsAbstractFile file, StsProgressPanel panel)
    {
        StsAsciiFile asciiFile = new StsAsciiFile(file);
        if(!asciiFile.openReadWithErrorMessage()) return false;

        StsCultureLine2D line = null;
        StsCulturePointSet2D pointSet = null;
        StsCultureText text = null;
        StringTokenizer st;
        String strToken;

        byte type = -1;
        int size = 1;
        int stroke = 0;

        try
        {
            StsColor color = StsColor.BLUE;
            StsColor overrideColor = cultureSetAttributes.panel.getStsColor();
            boolean override = cultureSetAttributes.panel.getOverrideColor();
            float coordinateScalar = cultureSetAttributes.panel.getScalar();        	
        	
            String aline = asciiFile.readLine();
            while(aline != null)
            {
                st = new StringTokenizer(aline, ", ");
                strToken = st.nextToken();
                if (strToken.equals(StsCultureObject2D.TYPES[StsCultureObject2D.LINE]))
                {
                    if(line != null) cultureSet2D.addObject(line); line = null;
                    if(text != null) cultureSet2D.addObject(text); text = null;
                    if(pointSet != null) cultureSet2D.addObject(pointSet); pointSet = null;
                    try
                    {
                        line = new StsCultureLine2D();
                        type = StsCultureObject2D.LINE;
                        if (st.countTokens() < 3)
                        {
                            panel.appendLine("Line type does not have 4 items, using defaults");
                            continue;
                        }
                        strToken = st.nextToken(); // Type of Line
                        strToken = st.nextToken(); // Width of Line
                        size = Integer.parseInt(strToken);
                        line.setSize(size);
                        strToken = st.nextToken();
                        color = StsColor.getColorByName(strToken);
                        if(override) color = overrideColor;
                        line.setStsColor(color); // Color of Line
                        strToken = st.nextToken();
                        stroke = Integer.parseInt(strToken);
                        line.setStroke(stroke);
                    }
                    catch (Exception e)
                    {
                        String fileLine = asciiFile.getLine();
                        StsException.systemError("Failed to read line for file line:\n" + fileLine);
                        panel.appendLine("ERROR: Failed to read line for file line:" + fileLine);
                        panel.setDescriptionAndLevel("Load failed", StsProgressBar.ERROR);
                        return false;
                    }
                }
                else if (strToken.equals(StsCultureObject2D.TYPES[StsCultureObject2D.TEXT]))
                {
                    if(line != null) cultureSet2D.addObject(line); line = null;
                    if(text != null) cultureSet2D.addObject(text); text = null;
                    if(pointSet != null) cultureSet2D.addObject(pointSet); pointSet = null;
                    try
                    {
                        text = new StsCultureText();
                        type = StsCultureObject2D.TEXT;
                        if (st.countTokens() < 3)
                        {
                            panel.appendLine("Text type does not have 4 items, using defaults");
                            continue;
                        }
                        strToken = st.nextToken(); // Type of Point
                        strToken = st.nextToken(); // Size of Point
                        size = Integer.parseInt(strToken);
                        if (size > 2)
                            text.setSize(size);
                        color = StsColor.getColorByName(st.nextToken());
                        if(override) color = overrideColor;
                        text.setStsColor(color); // Color of Point                        
                    }
                    catch (Exception e)
                    {
                        String fileLine = asciiFile.getLine();
                        StsException.systemError("Failed to read text set for line:\n" + fileLine);
                        panel.appendLine("ERROR: Failed to read line for file line:" + fileLine);
                        panel.setDescriptionAndLevel("Load failed", StsProgressBar.ERROR);
                        return false;
                    }
                }
                else if (strToken.equals(StsCultureObject2D.TYPES[StsCultureObject2D.POINT]))
                {
                    if(line != null) cultureSet2D.addObject(line); line = null;
                    if(text != null) cultureSet2D.addObject(text); text = null;
                    if(pointSet != null) cultureSet2D.addObject(pointSet); pointSet = null;
                    try
                    {
                        pointSet = new StsCulturePointSet2D();
                        type = StsCultureObject2D.POINT;
                        if (st.countTokens() < 3)
                        {
                            panel.appendLine("Pointset type does not have 4 items, using defaults\n");
                            continue;
                        }
                        strToken = st.nextToken(); // Type of Point
                        pointSet.setSymbolTypeByString(strToken);
                        strToken = st.nextToken(); // Size of Point
                        size = Integer.parseInt(strToken);
                        if (size > 2)
                            pointSet.setSize(size);
                        color = StsColor.getColorByName(st.nextToken());
                        if(override) color = overrideColor;
                        pointSet.setStsColor(color); // Color of Point
                    }
                    catch (Exception e)
                    {
                        String fileLine = asciiFile.getLine();
                        StsException.systemError("Failed to read point set for line:\n" + fileLine);
                        panel.appendLine("ERROR: Failed to read line for file line:" + fileLine);
                        panel.setDescriptionAndLevel("Load failed", StsProgressBar.ERROR);
                        return false;
                    }
                }
                else if (strToken.equals(StsCultureObject2D.TYPES[StsCultureObject2D.XY]))
                {
                    try
                    {
                        double x = Double.parseDouble(st.nextToken()) * coordinateScalar;
                        double y = Double.parseDouble(st.nextToken()) * coordinateScalar;
                        //float[] xy = model.getProject().getRotatedRelativeXY(x,y);

                        double z = model.getProject().getDepthMin();
                        double t = model.getProject().getTimeMin();

                        if (type == StsCultureObject2D.LINE)
                        {
                            if(st.hasMoreTokens())
                            {
                                if(model.getProject().getZDomain() == model.getProject().TD_DEPTH)
                                    z = Float.parseFloat(st.nextToken());
                                else
                                    t = Float.parseFloat(st.nextToken());
                            	line.setPlanar(false);                                
                            }
                            line.addPoint(x, y, z, t);
                        }
                        else if(type == StsCultureObject2D.POINT)
                        {
                            if(st.hasMoreTokens())
                            {
                                if(model.getProject().getZDomain() == model.getProject().TD_DEPTH)
                                    z = Double.parseDouble(st.nextToken());
                                else
                                    t = Double.parseDouble(st.nextToken());
                            	pointSet.setPlanar(false);
                            }
                            pointSet.addPoint(x, y, z, t);
                        }
                        else
                        {
                            if(st.countTokens() > 1) // Assume Z or T is supplied
                            {
                                if(model.getProject().getZDomain() == model.getProject().TD_DEPTH)
                                    z = Float.parseFloat(st.nextToken());
                                else
                                    t = Float.parseFloat(st.nextToken());
                                text.setPlanar(false);
                            }                        	
                            text.setPoint(x, y, z, t);
                            text.setText(st.nextToken());
                        }
                    }
                    catch (Exception e) {
                        String fileLine = asciiFile.getLine();
                        int nLines = asciiFile.getNLines();
                        StsException.systemError("Failed to parse float from line number " + nLines + ":\n" + fileLine);
                        panel.appendLine("Failed to parse float from line number " + nLines + ":\n" + fileLine + "\n");
                        panel.setDescriptionAndLevel("Load failed", StsProgressBar.ERROR);
                        return false;
                    }
                }
            aline = asciiFile.readLine();
            //System.out.println("Aline=" + aline);
            }
        }
        catch(Exception e)
        {
            String fileLine = asciiFile.getLine();
            panel.appendLine("ERROR: Failed to read line for file line:" + fileLine);
            panel.setDescriptionAndLevel("Load failed", StsProgressBar.ERROR);
            StsException.systemError("Failed to read line for file line:\n" + fileLine);
            return false;
        }

         // Add last point or line
         if(line != null)
         {
             cultureSet2D.addObject(line);
             line = null;
         }
         if(pointSet != null)
         {
             cultureSet2D.addObject(pointSet);
             pointSet = null;
         }
         if(text != null)
         {
             cultureSet2D.addObject(text);
             text = null;
         }
         if(cultureSet2D.getNumberOfObjects() == 0) 
         	return false;
         else
         	return true;
    }    
    public String getCultureSetName() { return cultureSetAttributes.panel.getName(); }
    public void setCultureSetName(String name)
    {
        cultureSetName = name;
    }
}
