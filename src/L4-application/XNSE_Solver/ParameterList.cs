using BoSSS.Foundation;
using Microsoft.SqlServer.Server;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver
{
    abstract class Parameter
    {
        public abstract IList<string> ParameterNames { get;}

        public DelParameterFactory Factory;

        public DelPartialParameterUpdate Update;
    }

    class ParameterList
    {
        List<Parameter> parameters;

        public ParameterList(int capacity = 10)
        {
            parameters = new List<Parameter>(capacity);
        }

        public void AddParameter(Parameter parameter)
        {
            parameters.Add(parameter);
        }

        public ICollection<DelParameterFactory> Factories(IList<string> names)
        {
            LinkedList<string> nameList = new LinkedList<string>(names);
            LinkedList<DelParameterFactory> parameterFactories = new LinkedList<DelParameterFactory>();
            //Find parameters and remove all found parameters from list;

            while(nameList.Count > 0)
            {
                string name = nameList.First.Value;
                nameList.RemoveFirst();
                //Find currentName
                for(int i = 0; i < parameters.Count; ++i)
                {
                    Parameter parameter = parameters[i];
                    if (parameter.ParameterNames.Contains(name))
                    {
                        if(parameter.Factory != null)
                        {
                            parameterFactories.AddLast(parameter.Factory);
                        }
                        foreach(string otherParamName in parameter.ParameterNames)
                        {
                            nameList.Remove(otherParamName);
                        }
                        break;
                    }
                }
            }
            return parameterFactories;
        }

        public ICollection<DelPartialParameterUpdate> ParameterUpdates(IList<string> names)
        {
            LinkedList<string> nameList = new LinkedList<string>(names);
            LinkedList<DelPartialParameterUpdate> parameterUpdates = new LinkedList<DelPartialParameterUpdate>();
            
            //Find parameters and remove all found parameters from list;
            while (nameList.Count > 0)
            {
                string name = nameList.First.Value;
                nameList.RemoveFirst();
                //Find currentName
                for (int i = 0; i < parameters.Count; ++i)
                {
                    Parameter parameter = parameters[i];
                    if (parameter.ParameterNames.Contains(name))
                    {
                        if(parameter.Update != null)
                        {
                            parameterUpdates.AddLast(parameter.Update);
                        }
                        foreach (string otherParamName in parameter.ParameterNames)
                        {
                            nameList.Remove(otherParamName);
                        }
                        break;
                    }
                }
            }
            return parameterUpdates;
        }
    }
}
