import React from "react";
import { Link } from "react-router-dom";
import { Tooltip, Box, Container, Typography, Card, CardContent, CardActions, Grid, IconButton } from '@mui/material';
import { GitHub, Description, OpenInNew, Info } from '@mui/icons-material';

const cardData = [
    {
        logo: 'cinemol.svg',
        title: 'CineMol',
        description: 'A web-based direct-to-SVG 3D small molecule drawer.',
        page: 'cinemol',
        links: {
            github: 'https://github.com/moltools/cinemol',
            doi: 'https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00851-y'
        }
    },
    {
        logo: 'paras.svg',
        title: 'PARAS(ECT)',
        description: 'Predict adenylation domain substrate specificity.',
        page: "https://paras.bioinformatics.nl/",
        links: {
            github: 'https://github.com/BTheDragonMaster/parasect',
            doi: 'https://www.biorxiv.org/content/10.1101/2025.01.08.631717v1'
        }
    },
    {
        logo: 'retromol.svg',
        title: 'RetroMol',
        description: 'Perform retrobiosynthetic analyses of modular natural products.',
        page: "https://retromol.bioinformatics.nl/",
        links: {
            github: 'https://github.com/moltools/RetroMol',
            // doi: ''
        }
    },
    {
        logo: 'biosynfoni.svg',
        title: 'Biosynfoni',
        description: 'Predict biosynthetic class for a given natural product.',
        page: "biosynfoni",
        links: {
            github: 'https://github.com/lucinamay/biosynfoni',
            // doi: ''
        }
    },
];

const LandingPage = () => {
    return (
        <Container 
            sx={{
                display: 'flex',
                flexDirection: 'column',
                padding: '2rem',
            }}
        >   
            <Box sx={{ display: 'flex', flexDirection: 'row', center: 'center', alignItems: 'center' }}>
                <img src="moltools.svg" alt="MolTools logo" style={{ width: '80px'}} />
                <Typography variant="h2" gutterBottom sx={{ margin: 0, paddingLeft: 2 }}>
                    MolTools (web app: v{process.env.REACT_APP_VERSION ? process.env.REACT_APP_VERSION : 'UNKNOWN'})
                </Typography>
            </Box>
            <Typography variant="body1" paragraph sx={{ paddingLeft: 2, marginBottom: '2rem' }}>
                Welcome to MolTools, a collection of tools for cheminformatics and molecular modeling of natural products.
            </Typography>
            <Grid container spacing={4}>
                {cardData.map((card, index) => (
                    <Grid item key={index} xs={12} sm={6}>
                        <Card
                            sx={{
                                ':hover': {
                                    boxShadow: 10, 
                                    transform: 'translateY(-5px)',
                                },
                                transition: '0.3s',
                                border: '1px solid #e0e0e0',
                                backgroundColor: '#fafafa',
                                boxShadow: 5,
                            }}
                            
                        >
                            <CardContent>
                                <Box sx={{ display: 'flex', flexDirection: 'row', alignItems: 'center', gap: '1rem' }}>
                                    <Box>
                                        <img src={card.logo} alt={card.title} style={{ width: '80px' }} />
                                    </Box>
                                    <Box>
                                        <Typography variant="h5" component="div">
                                            {card.title}
                                        </Typography>
                                        <Typography variant="body2" color="text.secondary">
                                            {card.description}
                                        </Typography>
                                    </Box>
                                </Box>
                            </CardContent>
                            <CardActions disableSpacing>
                                {card.page &&
                                    <Tooltip title="Launch application." arrow>
                                        <IconButton aria-label="open" component={Link} to={card.page}>
                                            <OpenInNew />
                                        </IconButton>
                                    </Tooltip>
                                }
                                {card.links.github &&
                                    <Tooltip title="View source code." arrow>
                                        <IconButton aria-label="github" component="a" href={card.links.github} target="_blank">
                                            <GitHub />
                                        </IconButton>
                                    </Tooltip>
                                }
                                {card.links.doi &&
                                    <Tooltip title="View publication." arrow>
                                        <IconButton aria-label="doi" component="a" href={card.links.doi} target="_blank">
                                            <Description />
                                        </IconButton>
                                    </Tooltip>
                                }
                            </CardActions>
                        </Card>
                    </Grid>
                ))}
            </Grid>
        </Container>
    );
};

export default LandingPage;