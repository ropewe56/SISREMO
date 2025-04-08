
# GET /public_power 
Public Power

Returns the public net electricity production for a given country for each production type. Subtype can be "solarlog" for Switzerland (ch).

Response schema:
```json
{
 "unix_seconds": list[int],
 "production_types": [
   {
       "name": str,
       "data": list[float]
   }
 ],
 "deprecated": bool
}
```

# GET /public_power_forecast
Public Power Forecast

Returns the forecast of the public net electricity production for a given country for each production type.

production_type: Can be solar, wind_onshore, wind_offshore or load
forecast_type: Can be current, intraday or day-ahead

If no dates are provided, values for today until forecast is available are returned. For load only the forecast type "day-ahead" is available.

Response schema:
```json
{
 "unix_seconds": list[int],
 "forecast_values": list[float],
 "production_type": str,
 "forecast_type": str,
 "deprecated": bool
}
```

# GET /total_power

Returns the total net electricity production (including industrial self supply) for a given country for each production type.

Currently only available for Germany.

Response schema:
```json
{
 "unix_seconds": list[int],
 "production_types": [
   {
       "name": str,
       "data": list[float]
   }
 ],
 "deprecated": bool
}
```

# GET installed_power

Returns the installed power for a specified country in GW except for battery storage capacity, which is given in GWh.

time_step: Time step can be either "yearly" or "monthly" (only for Germany)
installation_decommission: If true, the net installation / decommission numbers are returned instead of total installed power

Response schema:
```json
{
 "time": list[str],
 "production_types": [
   {
       "name": str,
       "data": list[float]
   }
 ],
 "deprecated": bool
}
```

# GET /cbet

Cross Border Electricity Trading

Returns the cross-border electricity trading (cbet) in GW between a specified country and its neighbors.

Positive values indicate an import of electricity, whereas negative values show electricity exports.

Response schema:
```json
{
    "unix_seconds": [int],
    "countries": [
        {
        "name": str,
        "data": [float]
        }
    ],
    "deprecated": bool
}
```

# GET /cbpf

Cross Border Physical Flows

Returns the cross-border physical flows (cbpfs) of electricity in GW between a specified country and its neighbors.

Positive values indicate an import of electricity, whereas negative values show electricity exports.

Response schema:
```json
{
    "unix_seconds": [int],
    "countries": [
        {
        "name": str,
        "data": [float]
        }
    ],
    "deprecated": bool
}
```