//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2008
//Author:       Larry K Wipperman
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Import;

import com.Sts.Actions.Wizards.OSWell.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.openspirit.*;
import com.openspirit.data.*;

import java.sql.*;
import java.util.*;

public class StsOpenSpiritImport
{
	private com.openspirit.OpenSpirit m_osp = null;
	private com.openspirit.user.SessionManager m_sessionMgr = null;
	private com.openspirit.user.Session m_session = null;
	private QueryScope m_queryScope = null;
	private QueryContext m_queryContext = null;
	private com.openspirit.metamodel.Model m_ospModel = null;
	private com.openspirit.unit.UnitCatalog m_unitCatalog = null;
	private com.openspirit.unit.UnitSystem m_unitSystem = null;
	private com.openspirit.carto.CoordinateReferenceSystem m_coordinateSystem = null;

    private StsModel model = null;
    private boolean debug = true;
	/**
	 * Default constructor ensures that an OpenSpirit installation exists on the
	 * server.
	 */
	public StsOpenSpiritImport(StsModel model)
	{
        this.model = model;
		try
		{
			// To run the example, an OpenSpirit runtime must be already
			// installed.
			String ospHome = com.openspirit.OpenSpiritFactory.getOpenSpiritHome().toString();
			if (ospHome == null)
			{
				outputMessage(StsMessage.ERROR, true,"StsOpenSpiritImport:The OSP_HOME Environment Variable has not been set.\n\n");
			}
		}
		catch (Exception e)
		{
			outputMessage(StsMessage.ERROR, true,"StsOpenSpiritImport:Failed to get the OSP_HOME Environment Variable:\nException: "
							+ e.getMessage() + "\n\n");
		}
	}

	/**
	 * initializeOpenSpirit creates an OpenSpirit instance, checks out an
	 * OpenSpirit Runtime license and probably an OspAppConnection license, and
	 * initializes some key OpenSpirit services.
	 */
	public boolean initializeOpenSpirit()
	{
		try
		{
			// Create an OpenSpirit instance
			m_osp = com.openspirit.OpenSpiritFactory.createOpenSpirit();
            if(debug) System.out.println("initializeOpenSpirit:Got instance...\n");
			// Check out the licenses. The exact content of the string parameter
			// doesn't matter. If another process uses the same string, another
			// OpenSpirit Runtime license will be checked out but no additional
			// OspAppConnection license will be checked out. If another process
			// has a different string, a new OspAppConnection license will be
			// checked out.
			m_osp.connect("StsApp");
            if(debug) System.out.println("initializeOpenSpirit:Connected to OpenSpirit instance...\n");
		}
		catch (InitializationException e)
		{
            e.printStackTrace();
			outputMessage(StsMessage.ERROR, true,"initializeOpenSpirit:InitializationException Connecting to OpenSpirit - Verify connection is available: "
							+ e.getMessage());
			return false;
		}
		catch (LicenseException e)
		{
			outputMessage(StsMessage.ERROR, true,"initializeOpenSpirit:LicenseException Connecting to OpenSpirit - Verify license is available: "
							+ e.getMessage());
			return false;
		}

		// The OpenSpirit model defines the data model to be used. We will use
		// the default
		// OpenSpirit data model
		try
		{
			m_ospModel = m_osp.getMetamodelService().getModel(com.openspirit.metamodel.MetamodelService.OPENSPIRIT_MODEL,null);
            if(debug) System.out.println("initializeOpenSpirit:Got OSP DataModel" + m_ospModel.getName());
		}
		catch (NotFoundException e)
		{
			outputMessage(StsMessage.ERROR, false,"initializeOpenSpirit:NotFoundException getting the Metamodel Service: " + e.getMessage() + "\n\n");
			exitOpenSpirit();
			return false;
		}
		catch (OspRuntimeException e)
		{
			outputMessage(StsMessage.ERROR, false,"initializeOpenSpirit:OspRuntimeException getting the Metamodel Service: "	+ e.getMessage() + "\n\n");
			exitOpenSpirit();
			return false;
		}

		// Initialize the unit catalog
		m_unitCatalog = m_ospModel.getUnitCatalog();

		return true;
	} // end of initializeOpenSpirit

	/**
	 * Get the com.openspirit.OpenSpirit instance. This is the basic connection
	 * to the OpenSpirit services.
	 * @return the com.openspirit.OpenSpirit instance.
	 */
	public com.openspirit.OpenSpirit getOpenSpiritInstance()
	{
		return m_osp;
	}

	/**
	 * Get the session names from the OpenSpirit session manager.
	 *
	 * @return a String array containing the session names. If no sessions
	 *         already exist the array will have zero length.
	 */
	public String[] getSessionNames()
	{
		// Initialize the session manager
		com.openspirit.user.UserPreferenceService userPrefService =	m_osp.getUserPreferenceService();
		com.openspirit.security.SecurityService securityService = m_osp.getSecurityService();
		m_sessionMgr = userPrefService.getSessionManager(securityService.getPrincipal());
        if(debug) System.out.println("getSessionNames:Got session manager " + m_sessionMgr.toString());

		String[] names = m_sessionMgr.getSessionNames();
        if(debug) System.out.println("getSessionNames:Found " + names.length + " sessions.");
		return names;
	}

	/**
	 * Set the session into the OpenSpirit instance and set the session's query
	 * scope and query context into the class variables.
	 *
	 * @param sessionName
	 *            is a String containing the name of the session to be set.
	 * @return true if successful or false if an error occurred.
	 */
	public boolean setSession(String sessionName)
	{
		if ((sessionName == null) || (sessionName.length() == 0))
		{
			outputMessage(StsMessage.ERROR, true,"setSession:The session name is null or empty.\n\n");
			return false;
		}

		try
		{
			m_session = m_sessionMgr.getSession(sessionName);
            if(debug) System.out.println("setSession:Successfully set session name to " + sessionName);
		}
		catch (NotFoundException e)
		{
			outputMessage(StsMessage.ERROR, false,"setSession:NotFoundException setting session "
					+ sessionName + ": " + e.getMessage() + "\n\n");
			return false;
		}

		// the query scope contains the projects set in the session. It is
		// possible that no projects are defined in the session.
		m_queryScope = m_session.getQueryScope();

		// the query context contains the Units and coordinate system set in the
		// session
		// Since we are just using well data, we will use Depth system. It is
		// possible that no units and/or coordinate system has been set.
		m_queryContext = m_session.getQueryContext(com.openspirit.user.VerticalDomain.DEPTH);

		return true;
	}

	/**
	 * Get the configured OpenSpirit data server names. These are normally the
	 * names given when the systems were configured for OpenSpirit. We will try
	 * to get them from the class queryScope first but it there are none, we
	 * will ask for them directly from the OpenSpirit DataService.
	 *
	 * @return a String array of data server names if they exist or an empty
	 *         array if they don't.
	 */
	public String[] getOspDataServers()
	{
		String[] servers = null;
		// first check if the class query scope has been set
        if(debug) System.out.println("getOspDataServers:Verify that the class query scope has been set.");
		if (m_queryScope != null && m_queryScope.getDataSources().length > 0)
		{
			servers = new String[m_queryScope.getDataSources().length];
            if(debug) System.out.println("getOspDataServers:Verified query scope and found " + servers.length + " servers.");
			DataSourceDefinition[] ds = m_queryScope.getDataSources();

			for (int i = 0; i < m_queryScope.getDataSources().length; i++)
			{
				servers[i] = ds[i].getName();
                //if(debug)
                    System.out.println("getOspDataServers:Server[" + i + "]= " + servers[i]);
			}
		}
		else
		{
			// get all of the data source types that support EpiWell_WellBore
			DataSourceTypeDefinition[] dst = null;
			try
			{
				dst = m_osp.getDataService().getDataSourceTypes(m_ospModel.getEntity("EpiWell_WellBore"));
			}
			catch (OspRuntimeException e)
			{
				outputMessage(StsMessage.ERROR, false,"getOspDataServers:OspRuntimeException getting OpenSpirit data source types: "
								+ e.getMessage() + "\n\n");
				return new String[0];
			}
			catch (NotFoundException e)
			{
				outputMessage(StsMessage.ERROR, false,"getOspDataServers:NotFoundException getting OpenSpirit data source types: "
								+ e.getMessage() + "\n\n");
				return new String[0];
			}

			// now we have to get the data servers for each type
			// we create strings containing the datastore type and version, followed
			// by the name.
			ArrayList<String> dservers = new ArrayList<String>();
             if(debug) System.out.println("getOspDataServers:Unable to verify query scope, retreived all data sources to support EpiWellbore");
			for (int i = 0; i < dst.length; i++)
			{
				String[] names = dst[i].getDataSourceNames();
				for (int j = 0; j < names.length; j++)
				{
					String typeName = dst[i].getName() + " " + dst[i].getVersion() +
										" - " + names[j];
                     if(debug) System.out.println("getOspDataServers:Adding data server named " + typeName + " to server list.");
					dservers.add(typeName);
				}
			}

			servers = new String[dservers.size()];
			if (dservers.size() > 0)
			{
				dservers.toArray(servers);
			}
		}

		return servers;
	} // end of getOspDataServers

	/**
	 * Get the project names for the input dataSource names.
	 *
	 * @param selectedDataSources
	 *            is a String array containing user selected data server names.
	 * @return a String array where each element is a String containing
	 *         dataServerName + " - " + projectName or an empty String array if
	 *         none are available or an error occurs.
	 */
	public String[] getProjectNames(String[] selectedDataSources)
	{
		String[] projectNames = null;
		int numServers = selectedDataSources.length;
        if(debug) System.out.println("getProjectNames:Retrieving project names from servers, Number of servers=" + numServers);
		if (numServers == 0)
		{
			outputMessage(StsMessage.ERROR, true,"getProjectNames:No selected OpenSpirit DataServers\n\n");
			return new String[0];
		}
		else
		{
			ArrayList<String> projNames = new ArrayList<String>();
			for (int i = 0; i < numServers; i++)
			{
				// parse the input string for the datastore name
				int iDash = selectedDataSources[i].indexOf("-");
				String name = selectedDataSources[i];
				if (iDash > -1)
					name = selectedDataSources[i].substring(iDash+2);
                if(debug)
                    System.out.println("getProjectNames:Getting data sources for server named=" + name);
				DataSourceDefinition[] dsd = null;
				try
				{
					dsd = m_osp.getDataService().getDataSources(name);
				}
				catch (NotFoundException e)
				{
					// don't worry about it, just keep going
					continue;
				}
				catch (OspRuntimeException e)
				{
					// this exception is more serious
					outputMessage(StsMessage.ERROR, true,"getProjectNames:OspRuntimeException getting OpenSpirit data sources for "
									+ "server name "
									+ selectedDataSources[i]
									+ ":\n" + e.getMessage() + "\n\n");
					return new String[0];
				}

				for (int j = 0; j < dsd.length; j++)
				{
					String[] names = dsd[j].getProjectNames();
                    if(debug) System.out.println("getProjectNames:" + names.length + " data sources found for server(" + name + ")");
					// if a datastore doesn't support projects, none will be returned.
					// However, the OpenSpirit datamodel will create "virtual" projects
					// if needed and we can query for them.
					if (names.length == 0)
					{
                        if(debug) System.out.println("getProjectNames:No projects found for server(" + name + ") we will query for virtual projects.");
						QueryScope queryScope;
						try
						{
							queryScope = m_osp.getDataService().getQueryFactory().createQueryScope(dsd[j]);
						}
						catch (BadArgumentsException e)
						{
							outputMessage(StsMessage.ERROR, false,"getProjectNames:BadArgumentsException creating queryScope:"
											+ e.getMessage() + "\n\n");
							return new String[0];
						}
						catch (OspRuntimeException e)
						{
							outputMessage(StsMessage.ERROR, false,"getProjectNames:OspRuntimeException creating queryScope:"
									+ e.getMessage() + "\n\n");
							return new String[0];
						}
				        OspConnection ospConnection = null;
				        OspPreparedStatement ospPreparedStatement = null;
				        OspResultSet ospResultSet = null;
				        try
				        {
				        	ospConnection = m_osp.getDataService().getOspConnection(queryScope, null);
							String query = "SELECT Name " +
					                       "FROM EpiProject_Project";
                            if(debug) System.out.println("getProjectNames:Query on data server " + name + " for available projects=" + query);
							ospPreparedStatement = ospConnection.prepareOspStatement(query, null);

							ospResultSet = ospPreparedStatement.executeOspQuery();

					        while (ospResultSet.next()) {
								String projName = ospResultSet.getString("Name");
								projNames.add(dsd[j].getName() + " - " + projName);
                                if(debug) System.out.println("getProjectNames:Adding project named " + projName + " from server named " + name);
							}

					        ospResultSet.close();
					        ospPreparedStatement.close();
					        ospConnection.close();
						}
				        catch (OspSQLException e)
				        {
				        	outputMessage(StsMessage.ERROR, true,"getProjectNames:OspSQLException processing OspResultSet:" + e.getMessage());
						}
				        catch (OspRuntimeException e)
				        {
				        	outputMessage(StsMessage.ERROR, false,"getProjectNames:OspRuntimeException:" + e.getMessage());
						}
				        catch (SQLException e)
				        {
				        	outputMessage(StsMessage.ERROR, false,"getProjectNames:SQLException processing OspResultSet:" + e.getMessage());
						}
						finally
						{
							try
							{
								if (ospResultSet != null)
									ospResultSet.close();
								if (ospPreparedStatement != null)
									ospPreparedStatement.close();
								if (ospConnection != null)
									ospConnection.close();
							}
							catch (SQLException e)
							{
								// do nothing
							}
						}
					}
					else
					{
                        if(debug) System.out.println("getProjectNames:" + names.length + " projects found for server(" + name + ")");
						for (int k = 0; k < names.length; k++)
						{
							projNames.add(dsd[j].getName() + " - " + names[k]);
                            if(debug) System.out.println("getProjectNames:Adding project (" + names[k] + ") to project list.");
						}
					}
				}
			}

			// extract the String array
			projectNames = new String[projNames.size()];
			projNames.toArray(projectNames);
		}

		return projectNames;
	} // end of getProjectNames

	/**
	 * Set the query scope to the array of input dataSource - projectName
	 * Strings.
	 *
	 * @param dataSourceProjectNames
	 *            is a String array of OpenSpirit DataSourceNames.
	 * @return true if successful or false if an error occurs.
	 */
	public boolean setQueryScope(String[] dataSourceProjectNames)
	{
		if (dataSourceProjectNames == null || dataSourceProjectNames.length == 0)
		{
			outputMessage(StsMessage.ERROR, true,"setQueryScope:No selected OpenSpirit Projects\n\n");
			return false;
		}
        if(debug) System.out.println("setQueryScope:Setting query scope to projects.");
		ArrayList<ProjectDefinition> projDefs =	new ArrayList<ProjectDefinition>(dataSourceProjectNames.length);
		DataSourceDefinition[] dsd = null;
		for (int i = 0; i < dataSourceProjectNames.length; i++)
		{
			// split the string into dataSourceName and projectName
			int iDash = dataSourceProjectNames[i].indexOf("-");
			String dataSourceName = dataSourceProjectNames[i].substring(0, iDash - 1);
			String projectName = dataSourceProjectNames[i].substring(iDash + 2);

            if(debug) System.out.println("setQueryScope:getting datasources for project " + projectName + " from datasoure " + dataSourceName);
			try
			{
				dsd = m_osp.getDataService().getDataSources(dataSourceName);

				// it is possible to be more than one data source with the same
				// name
				// for different data source types but I think it is pretty rare
				// in practice
				// so we will process however many come back
				for (int j = 0; j < dsd.length; j++)
				{
					projDefs.add(dsd[j].getProject(projectName));
				}
                 if(debug) System.out.println("setQueryScope:Adding project named " + projectName + " to project definitions.");
			}
			catch (NotFoundException e)
			{
				// don't worry about it, just keep going
				continue;
			}
			catch (OspRuntimeException e)
			{
				// this exception is more serious
				outputMessage(StsMessage.ERROR, true,"setQueryScope:OspRuntimeException getting OpenSpirit data sources for "
								+ "server name "
								+ dataSourceName
								+ ":\n"
								+ e.getMessage() + "\n\n");
				return false;
			}
		}

		ProjectDefinition[] projects = new ProjectDefinition[projDefs.size()];
		projDefs.toArray(projects);
        if(debug) System.out.println("setQueryScope:" + projDefs.size() + " projects defined within scope.");
		try
		{
			// it is possible that the projects array is empty if the datastore does not support projects
			// but in that case we can use the DataStoreDefinition for the query scope
			if (projects.length > 0)
            {
                            m_queryScope = m_osp.getDataService().
                                           getQueryFactory().createQueryScope(null,
                                projects);
                            if(debug) System.out.println("setQueryScope: Set query scope for project array:" + projects[0]);
            }
			else
            {
                            m_queryScope = m_osp.getDataService().
                                           getQueryFactory().createQueryScope(
                                dsd, null);
                            if(debug) System.out.println("setQueryScope: Set query scope for data service:" + dsd[0].getName());
            }
		}
		catch (BadArgumentsException e)
		{
			outputMessage(StsMessage.ERROR, false,"setQueryScope:BadArgumentsException creating a queryScope:\n"
							+ e.getMessage() + "\n\n");
			return false;
		}
		catch (OspRuntimeException e)
		{
			outputMessage(StsMessage.ERROR, false,"setQueryScope:OspRuntimeException creating a queryScope:\n"
							+ e.getMessage() + "\n\n");
			return false;
		}

		return true;
	} // end of setQueryScope

	/**
	 * Get the WellBoreLists from OSP with the existing queryScope, which must
	 * have already been set.
	 *
	 * @return a String[] of the WellBoreList names from the query scope. If
	 *         there are no lists to return or an error occurs, the String array
	 *         will be empty.
	 */
	public String[] getWellBoreLists()
	{
		OspConnection ospConnection = null;
		OspPreparedStatement ospPreparedStatement = null;
		OspResultSet ospResultSet = null;

		ArrayList<String> names = new ArrayList<String>();
		// open a data connection to OpenSpirit using the queryScope and
		// queryContext
		try
		{
			ospConnection = m_osp.getDataService().getOspConnection(m_queryScope, m_queryContext);

			String query = "SELECT Name FROM EpiWell_WellBoreList " +
						   "WHERE NumBores > 0";
            if(debug) System.out.println("getWellBoreLists:Retrieving wellbore list, query=" + query);
			ospPreparedStatement = ospConnection.prepareOspStatement(query,	null);

			// execute the query
            System.out.println("Submitting query for wellbore lists");
			ospResultSet = ospPreparedStatement.executeOspQuery();
            System.out.println("Returned from query for wellbore lists");
			while (ospResultSet.next())
			{
				names.add(ospResultSet.getString("Name"));
			}

			ospResultSet.close();
			ospPreparedStatement.close();
		}
		catch (OspSQLException e)
		{
			outputMessage(StsMessage.ERROR, true,"getWellBoreLists:OspSQLException getting Wellbore Lists:"
					+ e.getMessage() + "\n\n");
			return new String[0];
		}
		catch (OspRuntimeException e)
		{
			outputMessage(StsMessage.ERROR, false,"getWellBoreLists:OspRuntimeException getting Wellbore Lists:"
					+ e.getMessage() + "\n\n");
			return new String[0];
		}
		catch (SQLException e)
		{
			outputMessage(StsMessage.ERROR, false,"getWellBoreLists:SQLException getting Wellbore Lists:"
					+ e.getMessage() + "\n\n");
			return new String[0];
		} finally
		{
			// close open objects
			try
			{
				ospResultSet.close();
				ospPreparedStatement.close();
				ospConnection.close();
			}
			catch (SQLException e)
			{
				// ignore it
				return new String[0];
			}
		}

		String[] lists = new String[names.size()];
		names.toArray(lists);
        if(debug) System.out.println("getWellBoreLists:Found " + names.size() + " wellbores.");
		return lists;
	} // end of getWellBoreLists

	/**
	 * Get the well identifiers contained in the input wellbore lists
	 *
	 * @param wbLists
	 *            contains the user-selected wellbore list(s) to use.
	 * @return a String[] of "wellbore identifier - common name" strings for the
	 *         input wellbore list(s) or an empty string if an error occurs.
	 */
	public StsOSWell[] getWellsFromWellboreLists(String[] wbLists)
	{
            /*
		if (wbLists == null || wbLists.length == 0)
		{
			outputMessage(StsMessage.ERROR, true,"getWellsFromWellboreLists:The wellbore list array is empty\n\n");
			return new StsOSWell[0];
		}
*/
		OspConnection ospConnection = null;
		OspPreparedStatement ospPreparedStatement = null;
		OspResultSet ospResultSet = null;

		StsOSWell[] wells = null;
        if(debug) System.out.println("getWellsFromWellboreLists:Retriving wells from wellbore list...");
		try
		{
			ospConnection = m_osp.getDataService().getOspConnection(m_queryScope, m_queryContext);

			// this query doesn't work right now. There is a problem with the JOIN that I
			// am trying to get resolved with OpenSpirit. If worse comes to worse, I will
			// break this into two queries.
			String query = "SELECT Bores " +
						   "FROM EpiWell_WellBoreList " +
						   "WHERE Name IN (?)";

			ospPreparedStatement = ospConnection.prepareOspStatement(query,	null);
			ospPreparedStatement.setStringArray(1, wbLists);

			ospResultSet = ospPreparedStatement.executeOspQuery();
			ArrayList<DataKey[]> bores = new ArrayList<DataKey[]>();
			while (ospResultSet.next())
			{
				DataKey[] keys = ospResultSet.getDataKeyArray("Bores");
				bores.add(keys);
			}
			ospResultSet.close();
			ospPreparedStatement.close();

			// we want to put all the wellbore DataKeys into a single array
			// but first we need a total count of DataKeys
			int numBores = 0;
			for (int i=0; i<bores.size(); i++)
			{
				DataKey[] keys = bores.get(i);
				numBores += keys.length;
			}
            if(debug) System.out.println("getWellsFromWellboreLists:Wellbore count=" + numBores);

			DataKey[] allKeys = new DataKey[numBores];
			int index = 0;
			for (int i=0; i<bores.size(); i++)
			{
				DataKey[] keys = bores.get(i);
				for (int j=0; j<keys.length; j++)
					allKeys[index++] = keys[j];
			}
            if(debug) System.out.println("getWellsFromWellboreLists:Total keys=" + index);
			query = "SELECT wb.Identifier, wb.BoreAliases, wb.BoreAliasTypes, " +
						"w.Location, wb.PrimaryKey$ " +
			   		"FROM EpiWell_WellBore wb, EpiWell_Well w " +
			   		"WHERE wb.PrimaryKey$ IN (?)";

            if(debug) System.out.println("getWellsFromWellboreLists:Query to OpenSpirit=" + query);
			ospPreparedStatement = ospConnection.prepareOspStatement(query,	null);
			ospPreparedStatement.setDataKeyArray(1, allKeys);

			ospResultSet = ospPreparedStatement.executeOspQuery();
			wells = new StsOSWell[numBores];
			index = 0;
			while (ospResultSet.next())
			{
				String id = ospResultSet.getString("wb.Identifier");
				String[] aliases = ospResultSet.getStringArray("wb.BoreAliases");
				String[] types = ospResultSet.getStringArray("wb.BoreAliasTypes");

				// find the common wellbore name
				String commonName = null;
				for (int i = 0; i < aliases.length; i++)
				{
					if (types[i].toUpperCase().contains("COMMON"))
					{
						commonName = aliases[i];
						break;
					}
				}
                com.openspirit.data.spatial.Point surfaceLocation = (com.openspirit.data.spatial.Point)ospResultSet.getGeometry("w.Location");
                DataKey key = ospResultSet.getDataKey("wb.PrimaryKey$");
                StsOSWell well = new StsOSWell(id, commonName, surfaceLocation.getX(),	surfaceLocation.getY(), key);
                if(well != null)
               {
                   wells[index++] = well;
                   if(debug) System.out.println("getWellsFromWellboreLists:Added well #-" + index + " to wells array list" + well.toString());
               }
               else
                   System.out.println("getWellsFromWellboreLists:Failed to retrieve well (" + commonName + ") from wellbore in OpenSpirit.");
			}
			ospResultSet.close();
			ospPreparedStatement.close();
			ospConnection.close();
		}
		catch (OspSQLException e)
		{
			outputMessage(StsMessage.ERROR, true,"getWellsFromWellboreLists:OspSQLException getting wells for input wellbore lists:"
							+ e.getMessage() + "\n\n");
			return new StsOSWell[0];
		}
		catch (OspRuntimeException e)
		{
			outputMessage(StsMessage.ERROR, false,"getWellsFromWellboreLists:OspRuntimeException getting wells for input wellbore lists:"
							+ e.getMessage() + "\n\n");
			return new StsOSWell[0];
		}
		catch (SQLException e)
		{
			outputMessage(StsMessage.ERROR, false,"getWellsFromWellboreLists:SQLException getting wells for input wellbore lists:"
							+ e.getMessage() + "\n\n");
			return new StsOSWell[0];
		}
		finally
		{
			// close open objects
			try
			{
				if (ospResultSet != null)
					ospResultSet.close();
				if (ospPreparedStatement != null)
					ospPreparedStatement.close();
				if (ospConnection != null)
					ospConnection.close();
			}
			catch (SQLException e)
			{
				// ignore it
			}
		}

		return wells;
	} // end of getWellsFromWellboreLists

    /**
     * Get number of wells in project.
     *
     * @param projectName is the name of the project
     * @return int number of wells
     */
    public int getNumberOfWells(String projectName)
    {
        OspConnection ospConnection = null;
        OspPreparedStatement ospPreparedStatement = null;
        OspResultSet ospResultSet = null;

        int count = 0;
        try
        {
            ospConnection = m_osp.getDataService().getOspConnection(m_queryScope, m_queryContext);
            String query = "SELECT wb.Identifier FROM EpiWell_WellBore wb, EpiWell_Well w WHERE wb.Project$ = '" + projectName + "'";
            if (debug) System.out.println("getNumberOfWells:Executing query to retrieve wellbores, query=" + query);
            ospPreparedStatement = ospConnection.prepareOspStatement(query, null);
            ospResultSet = ospPreparedStatement.executeOspQuery();
            while (ospResultSet.next())
            {
                System.out.println("Total wells read = " + count++);
            }
        }
        catch (OspSQLException e)
        {
            outputMessage(StsMessage.ERROR, true,"getNumberOfWells:OspSQLException While getting wells for input wellbore lists.\n" +
                            "Verify that OpenSpirit licenses are available."	+ e.getMessage() + "\n\n");
            return 0;
        }
        catch (OspRuntimeException e)
        {
            outputMessage(StsMessage.ERROR, false,"getNumberOfWells:OspRuntimeException getting wells for input wellbore lists:"
                    + e.getMessage() + "\n\n");
            return 0;
        }
        catch (SQLException e)
        {
            outputMessage(StsMessage.ERROR, false,"getNumberOfWells:SQLException getting wells for input wellbore lists:"
                    + e.getMessage() + "\n\n");
            return 0;
        }
        return count;
    }
    /**
     * Get wells for the input project that are in the supplied list
     *
     * @param projectName is the name of the project
     * @param projectLimit is boolean used to limit list to project bounds
     * @param progressPanel is the progress panel used to display messages
     * @param wbList is the list of wellbore ids that we need to retrieve
     * @return an array of StsOSWell objects or an empty
     * array if no wells in the project or an error occurred.
     */
    public StsOSWell[] getProjectWellsInList(String projectName, boolean projectLimit, StsProgressPanel progressPanel, String[] wbList)
    {
        OspConnection ospConnection = null;
        OspPreparedStatement ospPreparedStatement = null;
        OspResultSet ospResultSet = null;

        int noSurfacePositionExcludedWells = 0, geographicalExcludedWells = 0, acceptedWells = 0;
        StsOSWell[] wells = null;

        ArrayList<StsOSWell> wbs = new ArrayList<StsOSWell>();
        try
        {
            for(int i=0; i<wbList.length; i++)
            {
                ospConnection = m_osp.getDataService().getOspConnection(m_queryScope, m_queryContext);
                String query =
                    "SELECT wb.BoreAliases, wb.BoreAliasTypes, " +
                    "w.Location, wb.PrimaryKey$ " +
                    "FROM EpiWell_WellBore wb, EpiWell_Well w " +
                    "WHERE wb.Identifier = '" + wbList[i] + "'";
                if (debug) System.out.println("getAllProjectWells:Executing query to retrieve wellbores, query=" +  query);
                ospPreparedStatement = ospConnection.prepareOspStatement(query, null);
                ospResultSet = ospPreparedStatement.executeOspQuery();
                ospResultSet.next();

                progressPanel.setDescription("Retrieved " + wbs.size() + " of " +
                                             (acceptedWells+noSurfacePositionExcludedWells+geographicalExcludedWells) +  " total wells.");
                if(debug)
                    System.out.println("getProjectWellsInList:Processing well #" + wbs.size());
                System.out.println(noSurfacePositionExcludedWells + " wells with no surface location, " +
                                   geographicalExcludedWells + " wells that are geographically excluded, and " +
                                   acceptedWells + " wells have been accepted.");

                //String id = ospResultSet.getString("wb.Identifier");
                String[] aliases = ospResultSet.getStringArray("wb.BoreAliases");
                String[] types = ospResultSet.getStringArray("wb.BoreAliasTypes");

                // find the common wellbore name
                String commonName = null;
                if(aliases != null && types != null)
                {
                    if (debug) System.out.println("getProjectWellsInList:Locating the wellbores common name from alias array of length " +
                                                  aliases.length);
                    for (int j = 0; j < aliases.length; j++)
                    {
                        if (types[j].toUpperCase().contains("COMMON"))
                        {
                            commonName = aliases[j];
                            break;
                        }
                    }
                }
                else
                {
                    commonName = "NoNameProvided";
                }
                if(debug) System.out.println("getProjectWellsInList:Common name is " + commonName);
                com.openspirit.data.spatial.Point surfaceLocation = (com.openspirit.data.spatial.Point)ospResultSet.getGeometry("w.Location");
                DataKey key = ospResultSet.getDataKey("wb.PrimaryKey$");
                if(debug) System.out.println("getProjectWellsInList:Datakey for well (" + commonName + ")=" + key.toString());
                StsOSWell well = null;
                if(surfaceLocation == null)
                {
                    noSurfacePositionExcludedWells++;
                    progressPanel.appendLine("Location is null for well (" + commonName + "), excluding from list.");
                    if(debug) System.out.println("getProjectWellsInList:SurfaceLocation is null for well (" + commonName + ") excluding from list.");
                    continue;
                }
                else
                {
                    if(debug) System.out.println("getProjectWellsInList:Creating StsOSWell (" + commonName + ")=" + key.toString());
                    if(model != null)
                    {
                        StsProject project = model.getProject();
                        if(projectLimit)
                        {
                            if(project.isInProjectBounds(surfaceLocation.getX(),surfaceLocation.getY()))
                            {
                                acceptedWells++;
                                well = new StsOSWell(wbList[i], commonName, surfaceLocation.getX(), surfaceLocation.getY(), key);
                                //System.out.println("getAllProjectWells:Well #" + wbs.size() + " named: "+ commonName + " has been created and added to list.");
                            }
                            else
                            {
                                geographicalExcludedWells++;
                                if(debug) System.out.println("getProjectWellsInList:Well " + commonName + " is outside project bounds, excluding from list");
                                continue;
                            }
                        }
                        else
                        {
                            acceptedWells++;
                            well = new StsOSWell(wbList[i], commonName, surfaceLocation.getX(), surfaceLocation.getY(), key);
                        }
                    }
                }
                if(well != null)
                {
                    wbs.add(well);
                    if(debug) System.out.println("getProjectWellsInList:Added well to wb#-" + wbs.size() + " list " + well.toString());
                }
                else
                    System.out.println("getProjectWellsInList:Failed to retrieve well (" + commonName + ") from OpenSpirit.");
            }
            progressPanel.appendLine("Retrieved " + wbs.size() + " wells out of a total of " + (acceptedWells+noSurfacePositionExcludedWells+geographicalExcludedWells) + " wells.");
            progressPanel.appendLine("    Rejected " + noSurfacePositionExcludedWells + " due to no surface location.");
            progressPanel.appendLine("    Rejected " + geographicalExcludedWells + " due to geographical exclusion.");
            if (ospResultSet != null)
                ospResultSet.close();
            if (ospPreparedStatement != null)
                ospPreparedStatement.close();
            if (ospConnection != null)
                ospConnection.close();
        }
        catch (OspSQLException e)
        {
            outputMessage(StsMessage.ERROR, true,"getProjectWellsInList:OspSQLException While getting wells for input wellbore lists.\n" +
                                    "Verify that OpenSpirit licenses are available."	+ e.getMessage() + "\n\n");
            return new StsOSWell[0];
        }
        catch (OspRuntimeException e)
        {
            outputMessage(StsMessage.ERROR, false,"getProjectWellsInList:OspRuntimeException getting wells for input wellbore lists:"
                            + e.getMessage() + "\n\n");
            return new StsOSWell[0];
        }
        catch (SQLException e)
        {
            outputMessage(StsMessage.ERROR, false,"getProjectWellsInList:SQLException getting wells for input wellbore lists:"
                            + e.getMessage() + "\n\n");
            return new StsOSWell[0];
        }
        finally
        {
            // close open objects
            try
            {
                if (ospResultSet != null)
                    ospResultSet.close();
                if (ospPreparedStatement != null)
                    ospPreparedStatement.close();
                if (ospConnection != null)
                    ospConnection.close();
            }
            catch (SQLException e)
            {
                // ignore it
            }
        }

        wells = new StsOSWell[wbs.size()];
        wbs.toArray(wells);
        System.out.println("getProjectWellsInList:Successfully retrieved " + wells.length + " wellbores from OpenSpirit.");
        return wells;
	}	// end of getProjectWellsInList
	/**
	 * Get all of the wells for the input project.
	 *
	 * @param projectName is the name of the project
	 * @return an array of StsOSWell objects or an empty
	 * array if no wells in the project or an error occurred.
	 */
	public StsOSWell[] getAllProjectWells(String projectName, boolean projectLimit, StsProgressPanel progressPanel)
	{
		OspConnection ospConnection = null;
		OspPreparedStatement ospPreparedStatement = null;
		OspResultSet ospResultSet = null;

        int noSurfacePositionExcludedWells = 0, geographicalExcludedWells = 0, acceptedWells = 0;
		StsOSWell[] wells = null;

		ArrayList<StsOSWell> wbs = new ArrayList<StsOSWell>();
		try
		{
			ospConnection = m_osp.getDataService().getOspConnection(m_queryScope, m_queryContext);
			String query = "SELECT wb.Identifier, wb.BoreAliases, wb.BoreAliasTypes, " +
								"w.Location, wb.PrimaryKey$ " +
							"FROM EpiWell_WellBore wb, EpiWell_Well w " +
                            "WHERE wb.Project$ = '" + projectName + "'";
            if(debug) System.out.println("getAllProjectWells:Executing query to retrieve wellbores, query=" + query);
			ospPreparedStatement = ospConnection.prepareOspStatement(query,	null);
			ospResultSet = ospPreparedStatement.executeOspQuery();
			while (ospResultSet.next())
			{
                progressPanel.setDescription("Retrieved " + wbs.size() + " of " +
                                                        (acceptedWells+noSurfacePositionExcludedWells+geographicalExcludedWells) +
                                                       " total wells.");

                if(debug)
                    System.out.println("getAllProjectWells:Processing well #" + wbs.size());
                System.out.println(noSurfacePositionExcludedWells + " wells with no surface location, " +
                                   geographicalExcludedWells + " wells that are geographically excluded, and " +
                                   acceptedWells + " wells have been accepted.");
				String id = ospResultSet.getString("wb.Identifier");
				String[] aliases = ospResultSet.getStringArray("wb.BoreAliases");
				String[] types = ospResultSet.getStringArray("wb.BoreAliasTypes");

				// find the common wellbore name
				String commonName = null;
                if(aliases != null && types != null)
                {
                    if (debug) System.out.println("getAllProjectWells:Locating the wellbores common name from alias array of length " +
                                                  aliases.length);
                    for (int i = 0; i < aliases.length; i++)
                    {
                        if (types[i].toUpperCase().contains("COMMON"))
                        {
                            commonName = aliases[i];
                            break;
                        }
                    }
                }
                else
                {
                    commonName = "NoNameProvided";
                }
                if(debug) System.out.println("getAllProjectWells:Common name is " + commonName);
                com.openspirit.data.spatial.Point surfaceLocation = (com.openspirit.data.spatial.Point)ospResultSet.getGeometry("w.Location");
                DataKey key = ospResultSet.getDataKey("wb.PrimaryKey$");
                if(debug) System.out.println("getAllProjectWells:Datakey for well (" + commonName + ")=" + key.toString());
                StsOSWell well = null;
                if(surfaceLocation == null)
                {
                    noSurfacePositionExcludedWells++;
                    progressPanel.appendLine("Location is null for well (" + commonName + "), excluding from list.");
                    if(debug) System.out.println("getAllProjectWells:SurfaceLocation is null for well (" + commonName + ") excluding from list.");
                    continue;
                }
                else
                {
                    if(debug) System.out.println("getAllProjectWells:Creating StsOSWell (" + commonName + ")=" + key.toString());
                    if(model != null)
                    {
                        StsProject project = model.getProject();
                        if(projectLimit)
                        {
                            if(project.isInProjectBounds(surfaceLocation.getX(),surfaceLocation.getY()))
                            {
                                acceptedWells++;
                                well = new StsOSWell(id, commonName, surfaceLocation.getX(), surfaceLocation.getY(), key);
                                //System.out.println("getAllProjectWells:Well #" + wbs.size() + " named: "+ commonName + " has been created and added to list.");
                            }
                            else
                            {
                                geographicalExcludedWells++;
                                if(debug) System.out.println("getAllProjectWells:Well " + commonName + " is outside project bounds, excluding from list");
                                continue;
                            }
                        }
                        else
                        {
                            acceptedWells++;
                            well = new StsOSWell(id, commonName, surfaceLocation.getX(), surfaceLocation.getY(), key);
                        }
                    }
                }
                if(well != null)
                {
                    wbs.add(well);
//                    progressPanel.setDescription("Retrieved " + wbs.size() + " of " +
//                                                 (acceptedWells+noSurfacePositionExcludedWells+geographicalExcludedWells) +
//                                                " total wells.");
                    if(debug) System.out.println("getAllProjectWells:Added well to wb#-" + wbs.size() + " list " + well.toString());
                }
                else
                    System.out.println("getAllProjectWells:Failed to retrieve well (" + commonName + ") from OpenSpirit.");
			}
            progressPanel.appendLine("Retrieved " + wbs.size() + " wells out of a total of " + (acceptedWells+noSurfacePositionExcludedWells+geographicalExcludedWells) + " wells.");
            progressPanel.appendLine("    Rejected " + noSurfacePositionExcludedWells + " due to no surface location.");
            progressPanel.appendLine("    Rejected " + geographicalExcludedWells + " due to geographical exclusion.");
			if (ospResultSet != null)
				ospResultSet.close();
			if (ospPreparedStatement != null)
				ospPreparedStatement.close();
			if (ospConnection != null)
				ospConnection.close();
		}
		catch (OspSQLException e)
		{
			outputMessage(StsMessage.ERROR, true,"getAllProjectWells:OspSQLException While getting wells for input wellbore lists.\n" +
                                    "Verify that OpenSpirit licenses are available."	+ e.getMessage() + "\n\n");
			return new StsOSWell[0];
		}
		catch (OspRuntimeException e)
		{
			outputMessage(StsMessage.ERROR, false,"getAllProjectWells:OspRuntimeException getting wells for input wellbore lists:"
							+ e.getMessage() + "\n\n");
			return new StsOSWell[0];
		}
		catch (SQLException e)
		{
			outputMessage(StsMessage.ERROR, false,"getAllProjectWells:SQLException getting wells for input wellbore lists:"
							+ e.getMessage() + "\n\n");
			return new StsOSWell[0];
		}
		finally
		{
			// close open objects
			try
			{
				if (ospResultSet != null)
					ospResultSet.close();
				if (ospPreparedStatement != null)
					ospPreparedStatement.close();
				if (ospConnection != null)
					ospConnection.close();
			}
			catch (SQLException e)
			{
				// ignore it
			}
		}

		wells = new StsOSWell[wbs.size()];
		wbs.toArray(wells);
        System.out.println("getAllProjectWells:Successfully retrieved " + wells.length + " wellbores from OpenSpirit.");
		return wells;
	}	// end of getAllProjectWells

	/**
	 * Create an StsOSWell object for each input wellbore data key.
	 * @return an StsWell array. The array is empty if an error occurs.
	 */
	public StsOSWell[] createStsOSWells(StsOSWell[] oswells, StsProgressPanel panel)
	{
        int cnt = 0;
		if (oswells == null || oswells.length == 0)
		{
			outputMessage(StsMessage.ERROR, true,"createStsOSWells:The input StsOSWell array is empty\n\n");
			return oswells;
		}
        if(debug) System.out.println("createStsOSWells:Creating " + oswells.length + " StsOsWells.");
		// put the wellbore data keys into an array
		DataKey[] keys = new DataKey[oswells.length];
		for (int i=0; i<oswells.length; i++)
			keys[i] = oswells[i].getDataKey();

		// we will create our own local query scope using the DataKeys
		QueryScope queryScope;
		try
		{
			queryScope = m_osp.getDataService().getQueryFactory().createQueryScope(keys);
		}
		catch (BadArgumentsException e)
		{
			outputMessage(StsMessage.ERROR, false,"createStsOSWells:BadArgumentsException creating queryScope:"
							+ e.getMessage() + "\n\n");
			return oswells;
		}
		catch (OspRuntimeException e)
		{
			outputMessage(StsMessage.ERROR, false,"createStsOSWells:OspRuntimeException creating queryScope:"
					+ e.getMessage() + "\n\n");
			return oswells;
		}

		OspConnection ospConnection = null;
		OspPreparedStatement ospPreparedStatement = null;
		OspResultSet ospResultSet = null;
//		StsWell[] stsWells = new StsWell[oswells.length];
		try
		{
            panel.appendLine("Establishing connection to selected datastore...");
			ospConnection = m_osp.getDataService().getOspConnection(queryScope, m_queryContext);

			// since we want to create an StsWell for each wellbore, we can loop through
			// each data key at a time.
			for (int i=0; i<keys.length; i++)
			{
                panel.appendLine("Accessing " + oswells[i].getName() + " information...");
				// get the well and wellbore information that we don't already have in the StsOSWell
				String query = "SELECT wb.IdentifierType, wb.Identifier, w.CurrentOperator, w.Field, " +
                                               "wb.spudDate, wb.CompletionDate, wb.Elevation, wb.ElevationDatum, " +
									"wb.PathMD, wb.PathTVD, wb.PathTWT, wb.PathXOffset, wb.PathYOffset, " +
									"wb.LogTraces, wb.Picks, wb.WellVelocities " +
								"FROM EpiWell_WellBore wb, EpiWell_Well w " +
								"WHERE wb.PrimaryKey$ = ?";
                if(debug) System.out.println("createStsOSWells:Query for key["+i+"]=" + query);

				ospPreparedStatement = ospConnection.prepareOspStatement(query,	m_queryContext);
				ospPreparedStatement.setDataKey(1, keys[i]);

				// execute the query
				ospResultSet = ospPreparedStatement.executeOspQuery();
                panel.appendLine("   Retrieved well information and trajectory.");
				String api = null;
				String uwi = null;
				String operator = null;
				String field = null;
				String date = null;
				float kbElev = 0.f;
 				float elev = 0.f;
                long spudDate = 0l;
                long completionDate = 0l;
                String elevDatum = null;
				float[] pathMD = null;
				float[] pathTVD = null;
				float[] pathTWT = null;
				float[] pathXOffset = null;
				float[] pathYOffset = null;
				DataKey[] logKeys = null;
				DataKey[] pickKeys = null;
				DataKey[] tdKeys = null;

				while (ospResultSet.next())
				{
					String idType = ospResultSet.getString("wb.IdentifierType");
					uwi = ospResultSet.getString("wb.Identifier");
					if (idType.toUpperCase().contains("API"))
					{
						api = ospResultSet.getString("wb.Identifier");
					}
					operator = ospResultSet.getString("w.CurrentOperator");
					field = ospResultSet.getString("w.Field");
                    Timestamp stamp = ospResultSet.getTimestamp("wb.CompletionDate");
                    if(stamp != null)
                    {
                        date = stamp.toString();
                        completionDate = stamp.getTime();
                    }
                    stamp = ospResultSet.getTimestamp("wb.SpudDate");
                    if(stamp != null)
                    {
                        spudDate = stamp.getTime();
                    }
                    if(debug) System.out.println("createStsOSWells:Retrieved well uwi=" + uwi + " Api=" + api);
					kbElev = ospResultSet.getFloat("wb.Elevation"); // Assumes elevation is from KB. Check elevation type.
                    elev = kbElev;
                    elevDatum = ospResultSet.getString("wb.ElevationDatum");
					com.openspirit.unit.FloatQuantitySeries pathMDSeries =
						ospResultSet.getFloatQuantitySeries("wb.PathMD");
					com.openspirit.unit.FloatQuantitySeries pathTVDSeries =
						ospResultSet.getFloatQuantitySeries("wb.PathTVD");
					com.openspirit.unit.FloatQuantitySeries pathTWTSeries =
						ospResultSet.getFloatQuantitySeries("wb.PathTWT");
					com.openspirit.unit.FloatQuantitySeries pathYOffsetSeries =
						ospResultSet.getFloatQuantitySeries("wb.PathYOffset");
					com.openspirit.unit.FloatQuantitySeries pathXOffsetSeries =
						ospResultSet.getFloatQuantitySeries("wb.PathXOffset");
					if (pathMDSeries != null)
					{
						pathMD = pathMDSeries.getValues();
						pathTVD = pathTVDSeries.getValues();
						// pathTWT is only available if there is a default well velocity
						if (pathTWT != null)
							pathTWT = pathTWTSeries.getValues();
						pathXOffset = pathXOffsetSeries.getValues();
						pathYOffset = pathYOffsetSeries.getValues();
					}
					logKeys = ospResultSet.getDataKeyArray("wb.LogTraces");
					pickKeys = ospResultSet.getDataKeyArray("wb.Picks");
					tdKeys = ospResultSet.getDataKeyArray("wb.WellVelocities");
				}
				if (ospResultSet != null)
					ospResultSet.close();
				if (ospPreparedStatement != null)
					ospPreparedStatement.close();

				/*********************************************
				 * add data to the current StsOSWell
				 *********************************************/
				oswells[i].setWellLabel(oswells[i].getName());
				oswells[i].setApi(api);
				oswells[i].setUwi(uwi);
				oswells[i].setOperator(operator);
				oswells[i].setField(field);
				oswells[i].setDate(date);
				oswells[i].setKbElev(kbElev);
                oswells[i].setElev(elev);
                oswells[i].setElevDatum(elevDatum);
                oswells[i].setSpudDate(spudDate);
                oswells[i].setCompletionDate(completionDate);
				if (pathMD != null)
				{
					StsFloatVector md = new StsFloatVector(pathMD);
					StsFloatVector tvd = new StsFloatVector(pathTVD);
					StsFloatVector xoffset = new StsFloatVector(pathXOffset);
					StsFloatVector yoffset = new StsFloatVector(pathYOffset);
					StsFloatVector twt = null;
					if (pathTWT != null)
						twt = new StsFloatVector(pathTWT);
					oswells[i].setDevData(md, tvd, twt, xoffset, yoffset);
				}

				if (pickKeys != null && pickKeys.length > 0)
				{
                    panel.appendLine("   Well has markers, retrieving...");
					// get the picks (markers) for this wellbore
					query = "SELECT PickName, PickType, Location, PickValueMD, " +
							"PickValueTVD, PickValueTWT " +
							"FROM EpiWell_WellPick " +
							"WHERE PrimaryKey$ IN (?) " +
	        				"ORDER BY PickValueMD";
                    if(debug) System.out.println("createStsOSWells:Query for well markers=" + query);
					ospPreparedStatement = ospConnection.prepareOspStatement(query,	m_queryContext);
					ospPreparedStatement.setDataKeyArray(1, pickKeys);

					ospResultSet = ospPreparedStatement.executeOspQuery();
                    cnt = 0;
					while (ospResultSet.next())
					{
						String pickName = ospResultSet.getString("PickName");
						String pickType = ospResultSet.getString("PickType");
                        if(debug) System.out.println("createStsOSWells:Adding well marker=" + pickName + " to well " + oswells[i].getName());
						byte type = 0;
						if (pickType.toUpperCase().contains("STRAT_MARKER"))
							type = 1;
						else if (pickType.toUpperCase().contains("FAULT_MARKER"))
							type = 3;
						com.openspirit.data.spatial.Point location =
							(com.openspirit.data.spatial.Point)ospResultSet.getGeometry("Location");
						float pickMd = ospResultSet.getFloat("PickValueMD");
						float pickTvd = ospResultSet.getFloat("PickValueTVD");
						// PickValueTWT will only have a value if there is a DefaultWellVelocity
						float pickTwt = ospResultSet.getFloat("PickValueTWT");

						//**************************************************
						// create new StsQSWellMarker and add to the StsOSWell
						//***************************************************
						StsOSWellMarker marker = new StsOSWellMarker(pickName, type,
								pickMd, pickTvd, pickTwt, location.getX(), location.getY());
						oswells[i].addWellMarker(marker);
                        cnt++;
					}
					if (ospResultSet != null)
						ospResultSet.close();
					if (ospPreparedStatement != null)
						ospPreparedStatement.close();
                    panel.appendLine("   Retrieved " + cnt +" well markers from datastore.");
				}

				if (logKeys != null && logKeys.length > 0)
				{
                    panel.appendLine("   Well has logs, retrieving...");
					// get the logs for this wellbore
					query = "SELECT TraceName, TraceKind, NullValue, TraceType, TraceVersion, " +
								"MinTraceValue, MaxTraceValue, " +
								"TraceData, TraceIndexMD, TraceIndexTVD " +
							"FROM EpiWell_WellLogTrace " +
							"WHERE PrimaryKey$ IN (?)";
                    if(debug) System.out.println("createStsOSWells:Query for well logs=" + query);
					ospPreparedStatement = ospConnection.prepareOspStatement(query,	m_queryContext);
					ospPreparedStatement.setDataKeyArray(1, logKeys);

					ospResultSet = ospPreparedStatement.executeOspQuery();
                    cnt = 0;
					while (ospResultSet.next())
					{
						String name = ospResultSet.getString("TraceName");
						String kind = ospResultSet.getString("TraceKind");
						float nullValue = ospResultSet.getFloat("NullValue");
						String sampType = ospResultSet.getString("TraceType");
						int version = ospResultSet.getInt("TraceVersion");
						float minValue = ospResultSet.getFloat("MinTraceValue");
						float maxValue = ospResultSet.getFloat("MaxTraceValue");
		        		com.openspirit.unit.FloatQuantitySeries traceDataSeries =
							ospResultSet.getFloatQuantitySeries("TraceData");
		        		com.openspirit.unit.FloatQuantitySeries indexMdSeries =
							ospResultSet.getFloatQuantitySeries("TraceIndexMD");
		        		com.openspirit.unit.FloatQuantitySeries indexTvdSeries =
							ospResultSet.getFloatQuantitySeries("TraceIndexTVD");
                        if(debug) System.out.println("createStsOSWells:Found log with name=" + name);
						float[] logValues = null;
						float[] md = null;
						float[] tvd = null;
						if (traceDataSeries != null)
							logValues = traceDataSeries.getValues();
						if (indexMdSeries != null)
							md = indexMdSeries.getValues();
						if (indexTvdSeries != null)
							tvd = indexTvdSeries.getValues();

						//**************************************************
						// create a new StsLogCurve and add to the StsOSWell
						//***************************************************

                        //StsLogVector mdVector = new StsLogVector((byte)3, md);
                        //StsLogVector tvdVector = new StsLogVector((byte)2, tvd);

                        StsLogVector[] vector = StsWellKeywordIO.constructLogVectors(new String[] {StsLogVector.types[StsLogVector.MDEPTH]},
                                                               oswells[i].getName(), StsLogVector.WELL_LOG_PREFIX);
                        StsLogVector mdVector = vector[0];
                        mdVector.setValues(new StsFloatVector(md));

                        vector = StsWellKeywordIO.constructLogVectors(new String[] {StsLogVector.types[StsLogVector.DEPTH]},
                                                                             oswells[i].getName(), StsLogVector.WELL_LOG_PREFIX);
						StsLogVector tvdVector = vector[0];
                        if(tvdVector != null)
                        {
                            if(tvd != null)
                                tvdVector.setValues(new StsFloatVector(tvd));
                            else
                                // No tvd curve so use MD ---- Not sure this is valid but crashes without it. Should construct vertical well path
                                tvdVector.setValues(new StsFloatVector(md));
                        }
						vector = StsWellKeywordIO.constructLogVectors(new String[] {name}, oswells[i].getName(), StsLogVector.WELL_LOG_PREFIX);
                        vector[0].setNullValue(nullValue); // Must set prior to setValues to reset user null to S2S null.
                        vector[0].setValues(new StsFloatVector(logValues));
                        vector[0].setMinMaxAndNulls(nullValue);

                        if((mdVector.getValues().getSize() != tvdVector.getValues().getSize()) && (mdVector.getValues().getSize() != vector[0].getValues().getSize()))
                            System.out.println("Number of values does not match between log, tvd and md...\n");

						StsLogCurve log = new StsLogCurve(mdVector, tvdVector, vector[0], version);
						oswells[i].addLog(log);
                        if(debug) System.out.println("createStsOSWells:Added log (" + name + ") to well " + oswells[i].getName());
                        cnt++;
					}
					if (ospResultSet != null)
						ospResultSet.close();
					if (ospPreparedStatement != null)
						ospPreparedStatement.close();
                    panel.appendLine("   Retrieved " + cnt +" logs from datastore.");
				}

                cnt = 0;
				if (tdKeys != null && tdKeys.length > 0)
				{
                    panel.appendLine("   Well has time-depth curves, retrieving...");
					// if we have a default well velocity, just get it.
					// otherwise, get all available velocities
					query = "SELECT WellVelocityName, Datum, DisplayShift, SeismicReferenceDatum, " +
								"VelocityToSeismicDatum, TimeDepthMD, TimeDepthTVD, TimeDepthTWT, " +
								"IsBoreDefault, LastModifiedDate " +
							"FROM EpiWell_WellVelocity " +
							"WHERE PrimaryKey$ IN (?)";
                    if(debug) System.out.println("createStsOSWells:Query for well td curves=" + query);
					ospPreparedStatement = ospConnection.prepareOspStatement(query,	m_queryContext);
					ospPreparedStatement.setDataKeyArray(1, tdKeys);

					ospResultSet = ospPreparedStatement.executeOspQuery();
					java.sql.Timestamp lastTime = null;
					lastTime.setTime(0);
					while (ospResultSet.next())
					{
						String name = ospResultSet.getString("WellVelocityName");
						float datum = ospResultSet.getFloat("Datum");	// depth where time=0
						float displayShift = ospResultSet.getFloat("DisplayShift");	// static shift for display
						float seismicDatum = ospResultSet.getFloat("SeismicReferenceDatum");	// seismic reference elevation
						float srdVelocity = ospResultSet.getFloat("VelocityToSeismicDatum");
		        		com.openspirit.unit.FloatQuantitySeries tdMdSeries =
							ospResultSet.getFloatQuantitySeries("TimeDepthMD");
		        		com.openspirit.unit.FloatQuantitySeries tdTvdSeries =
							ospResultSet.getFloatQuantitySeries("TimeDepthTVD");
		        		com.openspirit.unit.FloatQuantitySeries tdTimeSeries =
							ospResultSet.getFloatQuantitySeries("TimeDepthTWT");
		        		float[] tdMd = null;
		        		float[] tdTvd = null;
		        		float[] tdTime = null;
		        		if (tdMdSeries != null)
		        		{
		        			tdMd = tdMdSeries.getValues();
		        			tdTvd = tdTvdSeries.getValues();
		        			tdTime = tdTimeSeries.getValues();
		        		}
		        		boolean isDefault = ospResultSet.getBoolean("IsBoreDefault");
		        		java.sql.Timestamp time = ospResultSet.getTimestamp("LastModifiedDate");

		        		// if there is more than one well velocity, we just want either
		        		// the default or the newest
		        		if (!isDefault && time.before(lastTime))
                        {
                            if(debug) System.out.println("createStsOSWells:More than one td curve was found, skipping current one to loaded latest.");
                            continue;
                        }
						//**************************************************
						// create a new StsLogCurve and add to the StsOSWell
						//***************************************************
                        StsLogVector[] vector = StsWellKeywordIO.constructLogVectors(new String[] {StsLogVector.types[StsLogVector.MDEPTH]},
                            oswells[i].getName(), StsLogVector.WELL_TD_PREFIX);
						StsLogVector mdVector = vector[0];
                        vector[0].setValues(new StsFloatVector(tdMd));

                        vector = StsWellKeywordIO.constructLogVectors(new String[] {StsLogVector.types[StsLogVector.DEPTH]},
                            oswells[i].getName(), StsLogVector.WELL_TD_PREFIX);
						StsLogVector tvdVector = vector[0];
                        vector[0].setValues(new StsFloatVector(tdTvd));

                        vector = StsWellKeywordIO.constructLogVectors(new String[] {StsLogVector.types[StsLogVector.TWT]},
                            oswells[i].getName(), StsLogVector.WELL_TD_PREFIX);
                        vector[0].setValues(new StsFloatVector(tdTime));

						StsLogCurve tdCurve = new StsLogCurve(mdVector, tvdVector, vector[0], 0);
						oswells[i].setTdCurve(tdCurve);
                        if(debug) System.out.println("createStsOSWells:Added tdCurve to well " + oswells[i].getName());
						lastTime = time;
                        cnt++;
					}
					if (ospResultSet != null)
						ospResultSet.close();
					if (ospPreparedStatement != null)
						ospPreparedStatement.close();
                    panel.appendLine("   Retrieved " + cnt +" time-depth curves from datastore. Using final curve.");
				}

			}	// end of foreach wellbore DataKey
		}
		catch (OspSQLException e)
		{
			outputMessage(StsMessage.ERROR, true,"createStsOSWells:OspSQLException getting Wellbore data:"
					+ e.getMessage() + "\n\n");
		}
		catch (OspRuntimeException e)
		{
			outputMessage(StsMessage.ERROR, false,"createStsOSWells:OspRuntimeException getting Wellbore data:"
							+ e.getMessage() + "\n\n");
		}
		catch (SQLException e)
		{
			outputMessage(StsMessage.ERROR, false,"createStsOSWells:SQLException getting Wellbore data:"
					+ e.getMessage() + "\n\n");
		}
		finally
		{
			// close open objects
			try
			{
				if (ospResultSet != null)
					ospResultSet.close();
				if (ospPreparedStatement != null)
					ospPreparedStatement.close();
				if (ospConnection != null)
					ospConnection.close();
			}
			catch (SQLException e)
			{
				// ignore it
			}
		}
		return oswells;
	} // end of createStsWells from StsOSWells

	/**
	 * exitOpenSpirit disconnects the OpenSpirit instance and checks the
	 * OspAppConnection and OpenSpirit Runtime licenses back in.
	 */
	public void exitOpenSpirit()
	{
		if ((m_osp != null) && m_osp.isConnected())
		{
			m_osp.disconnect();
		}
	}

	public void outputMessage(int level, boolean dialog, String msg)
	{
		if(dialog)
			new StsMessage(model.win3d, level, msg);
		else
			StsMessage.printMessage(msg);
		StsMessageFiles.logMessage(msg);
	}
}
