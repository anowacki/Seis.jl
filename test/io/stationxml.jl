using Test
using Dates: DateTime
using Seis
import Seis.SeisStationXML
import StationXML

# StationXML IO testing is done by the StationXML.jl package;
# these tests are for the Seis interface only.
@testset "StationXML" begin
    @testset "SeisStationXML" begin
        @testset "filter_xml" begin
            sxml = StationXML.FDSNStationXML(
                source="Me", created=DateTime("2000-01-01"), schema_version="1.1"
            )
            net1 = StationXML.Network(code="AN")
            net2 = StationXML.Network(code="BO")
            sta1 = StationXML.Station(
                code="FAKE", longitude=1, latitude=2, elevation=3,
                site=StationXML.Site(name="Fake site")
            )
            sta2 = StationXML.Station(
                code="FAKE2", longitude=2, latitude=3, elevation=4,
                site=StationXML.Site(name="Fake site 2")
            )
            cha1 = StationXML.Channel(
                code="UHT", longitude=10, latitude=20, elevation=30,
                depth=40, location_code=""
            )
            cha2 = StationXML.Channel(
                code="XYZ", longitude=20, latitude=30, elevation=40,
                depth=50, location_code="00"
            )

            @testset "Channel" begin
                sxml_all_channels = deepcopy(sxml)
                sxml_one_channel = deepcopy(sxml)
                append!(sxml_all_channels.network, deepcopy([net1, net2]))
                for n in sxml_all_channels.network
                    append!(n.station, deepcopy([sta1, sta2]))
                    for s in n.station
                        append!(s.channel, deepcopy([cha1, cha2]))
                    end
                end
                push!(sxml_one_channel.network, deepcopy(net1))
                push!(sxml_one_channel.network[1].station, deepcopy(sta1))
                push!(sxml_one_channel.network[1].station[1].channel, deepcopy(cha1))
                @test SeisStationXML.filter_stationxml(
                    sxml_all_channels, net1, sta1, cha1
                ) == sxml_one_channel

                @test sxml_one_channel.network[1].selected_num_stations === missing
                @test sxml_one_channel.network[1].station[1].selected_num_channels === missing
            end

            @testset "Selected number channels/stations" begin
                net = StationXML.Network(
                    code="AN", selected_number_stations=2, station=[sta1, sta2]
                )
                sxml′ = StationXML.FDSNStationXML(
                    source="Me", created=DateTime("2000-01-01"), schema_version="1.1",
                    network=[net]
                )
                sxml_filtered = SeisStationXML.filter_stationxml(sxml′, net, sta1, cha1)
                @test sxml_filtered.network[1].selected_number_stations == 1
            end
        end


        @testset "_getifnotmissing" begin
            @test SeisStationXML._getifnotmissing((a=1, b=missing), :a) == 1
            @test SeisStationXML._getifnotmissing(missing, :a) === missing
        end
    end

    @testset "Reading" begin
        file = joinpath(dirname(pathof(StationXML)), "..", "test", "data", "JSA.xml")
        
        @testset "File" begin
            @testset "Station{$T}" for T in (Float32, Float64)
                @test read_stationxml(
                    file, GeogStation{T}; warn=false, full=false
                ) isa Vector{GeogStation{T}}
            end

            @testset "Default type" begin
                @test read_stationxml(file) isa Vector{GeogStation{Float64}}
            end

            @testset "Properties" begin
                stas = read_stationxml(file)
                @test length(stas) == 6
                @test stas.net == fill("GB", 6)
                @test stas.sta == fill("JSA", 6)
                @test stas.loc == fill("", 6)
                @test stas.cha == ["BHE", "BHN", "BHZ", "HHE", "HHN", "HHZ"]

                sta = stas[1]
                @test sta.lon == -2.171698
                @test sta.lat == 49.187801
                @test sta.elev == 39.0
                @test sta.azi == 90.0
                @test sta.inc == 90.0
                @test sta.meta.startdate == DateTime("2007-09-06T00:00:00")
                @test sta.meta.enddate == DateTime("2599-12-31T23:59:59")
                @test sta.meta.burial_depth == 0.0
            end
        end

        @testset "IO" begin
            @test_throws(
                "XMLError: Document is empty from XML parser (code: 4, line: 1)",
                read_stationxml(IOBuffer(""))
            )
            @test read_stationxml(IOBuffer(read(file, String))) == read_stationxml(file)
        end

        @testset "String" begin
            @test_throws(ArgumentError, parse_stationxml(""))
            @test_throws(
                "XMLError: Premature end of data in tag FDSNStationXML line 2 from XML parser (code: 77, line: 7)",
                parse_stationxml(
                    """
                    <?xml version="1.0" encoding="UTF-8"?>
                    <FDSNStationXML xmlns="http://www.fdsn.org/xml/station/1" schemaVersion="1.0">
                      <Source>StationXML.jl</Source>
                      <Sender>Test</Sender>
                      <Created>2019-04-01T09:58:21.123456789+01:00</Created>
                    <!-- </FDSNStationXML> -->
                    """
                )
            )
            @test parse_stationxml(read(file, String)) == read_stationxml(file)
        end
    end
end